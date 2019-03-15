from headers import *

##################################################################################
##################################################################################

class Catalog(object):

   def __init__(self, U, MassConversion, name="test", nameLong=None, pathInCatalog="", save=False):

      self.U = U
      self.MassConversion = MassConversion
      self.name = name
      if nameLong is None:
         self.nameLong = self.name
      else:
         self.nameLong = nameLong
      self.pathInCatalog = pathInCatalog
      
      # Output path
      self.pathOut = "./output/catalog/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      # catalog path
      self.pathOutCatalog = self.pathOut + "/catalog.txt"
      
      # Figures path
      self.pathFig = "./figures/catalog/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      
      if save:
         self.readInputCatalog()
         self.addHaloMass()
         self.addIntegratedTau()
         self.addIntegratedKSZ()
         self.addIntegratedY()
         self.writeCatalog()

      self.loadCatalog()
   

   ##################################################################################
   ##################################################################################

   def copy(self, name="test", nameLong=None):
      """Copy a catalog class, with the option of changing the name.
      """
      # First copy the output catalog
      # new catalog path
      newPathOut = "./output/catalog/"+name
      if not os.path.exists(newPathOut):
         os.makedirs(newPathOut)
      newPathOutCatalog = newPathOut + "/catalog.txt"
      # copy the output catalog
      copyfile(self.pathOutCatalog, newPathOutCatalog)
      
      # Then copy the catalog properties
      newCat = Catalog(self.U, self.MassConversion, name=name, nameLong=nameLong, pathInCatalog=self.pathInCatalog, save=False)
      return newCat


   ##################################################################################
   ##################################################################################

   def readInputCatalog(self):
      print "- read input catalog from "+self.pathInCatalog
      data = np.genfromtxt(self.pathInCatalog)
      self.nObj = len(data[:,0])
      #
      # sky coordinates and redshift
      self.RA = data[:,0] # [deg]
      self.DEC = data[:,1]   # [deg]
      self.Z = data[:,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:,3]   # [Mpc/h]
      self.coordY = data[:,4]   # [Mpc/h]
      self.coordZ = data[:,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:,6]   # [Mpc/h]
      self.dY = data[:,7]   # [Mpc/h]
      self.dZ = data[:,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:,10]   # [Mpc/h]
      self.dZKaiser = data[:,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:,12]   #[km/s]
      self.vY = data[:,13]   #[km/s]
      self.vZ = data[:,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:,15]  # [km/s]   from spherical catalo
      self.vTheta = data[:,16]   # [km/s]
      self.vPhi = data[:,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:,18]   # [M_sun], from Maraston et al


   ##################################################################################

   def addHaloMass(self):
      """Generate halo masses in M_sun,
      from stellar masses in M_sun.
      """
      print "- add halo masses"
      # flag: 1 if object has mass
      self.hasM = np.zeros(self.nObj)
      
      self.Mvir = np.zeros(self.nObj)
      for iObj in range(self.nObj):
         mStellar = self.Mstellar[iObj]
         if (mStellar>1.e3) and not np.isnan(mStellar):
            self.hasM[iObj] = True
            self.Mvir[iObj] = self.MassConversion.fmStarTomVir(mStellar)

      # for object without a mass, use the mean mass from the others
      if np.sum(self.hasM)>0:
         meanMstellar = np.mean(self.Mstellar[self.hasM.astype('bool')])
         self.Mstellar[~self.hasM.astype('bool')] = meanMstellar
         #
         meanMvir = np.mean(self.Mvir[self.hasM.astype('bool')])
         self.Mvir[~self.hasM.astype('bool')] = meanMvir
      # if no object has a mass, make a random guess, rather than keeping 0
      else:
         self.Mstellar = 2.6e11   # random guess
         self.Mvir = self.MassConversion.fmStarTomVir(2.6e11)



   def addIntegratedTau(self):
      """integrated optical depth to Thompson scattering: int d^2theta n_e^2d sigma_T
      = (total nb of electrons) * sigma_T / (a chi)^2
      [sr]
      """
      print "- add integrated tau"
      # convert from total mass to baryon mass
      # assuming cosmological abundance of baryons, and all baryons are in gas
      self.integratedTau = self.Mvir * self.U.bg.Omega0_b/self.U.bg.Omega0_m
      
      # convert from total baryon mass to electron total number
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      msun = 1.989e30   # solar mass (kg)
      factor = (me + nH_ne*mH + nHe_ne*mHe) * (1./msun)   # total mass per electron in (Msun)
      self.integratedTau /= factor
      
      # multiply by Thomson cross section (physical)
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      self.integratedTau *= sigma_T / (mpc / self.U.bg.h)**2
      
      # divide by (a chi)^2
      self.integratedTau /= (self.U.bg.comoving_distance(self.Z) / (1.+self.Z))**2

   
   def addIntegratedKSZ(self):
      """Integrated kSZ signal: int d^2theta n_e sigma_T v/c Tcmb
      in muK * sr
      """
      print "- add integrated kSZ"
      self.integratedKSZ = - self.integratedTau * (self.vR/3.e5) * 2.726e6

   
   def addIntegratedY(self, nu=150.e9):
      """Integrated tSZ signal: int d^2theta n_e sigma_T (k_B T_e / m_e c^2)
      in sr.
      To get dT in muK*sr, multiply by Tcmb * f(nu).
      Simple power-law fit to Greco et al 2014, fig4.
      """
      print "- add integrated y"
      
      # in arcmin^2
      yCcyltilda = (self.Mstellar/1.e11)**3.2 * 1.e-6

      # in arcmin^2
      yCcyl = yCcyltilda * (self.U.hubble(self.Z) / self.U.hubble(0.))**(2./3.)
      yCcyl /= (self.U.bg.comoving_distance(self.Z) / (500.*self.U.bg.h))**2
      # in sr
      yCcyl *= (np.pi/180./60.)**2
      
      self.integratedY = yCcyl


   ##################################################################################

   def writeCatalog(self):
      print "- write full catalog to "+self.pathOutCatalog
      data = np.zeros((self.nObj,24))
      #
      # sky coordinates and redshift
      data[:,0] = self.RA # [deg]
      data[:,1] = self.DEC   # [deg]
      data[:,2] = self.Z
      #
      # observed cartesian coordinates
      data[:,3] = self.coordX   # [Mpc/h]
      data[:,4] = self.coordY   # [Mpc/h]
      data[:,5] = self.coordZ   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      data[:,6] = self.dX   # [Mpc/h]
      data[:,7] = self.dY   # [Mpc/h]
      data[:,8] = self.dZ   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      data[:,9] = self.dXKaiser   # [Mpc/h] from cartesian catalog difference
      data[:,10] = self.dYKaiser   # [Mpc/h]
      data[:,11] = self.dZKaiser   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      data[:,12] = self.vX   #[km/s]
      data[:,13] = self.vY   #[km/s]
      data[:,14] = self.vZ   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      data[:,15] = self.vR  # [km/s]   from spherical catalo
      data[:,16] = self.vTheta   # [km/s]
      data[:,17] = self.vPhi  # [km/s]
      #
      # Stellar mass
      data[:,18] = self.Mstellar   # [M_sun], from Maraston et al
      #
      # Halo mass
      data[:,19] = self.hasM  # flag=1 if mass is known
      data[:,20] = self.Mvir   # [M_sun]
      #
      # Integrated optical depth [sr]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      data[:,21] = self.integratedTau   # [sr]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e^2d sigma_T v/c Tcmb
      data[:, 22] = self.integratedKSZ # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e^2d sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      data[:, 23] = self.integratedY # [sr]
      #
      np.savetxt(self.pathOutCatalog, data)


   def loadCatalog(self):
      print "- load full catalog from "+self.pathOutCatalog
      data = np.genfromtxt(self.pathOutCatalog)
      self.nObj = len(data[:,0])
      #
      # sky coordinates and redshift
      self.RA = data[:,0] # [deg]
      self.DEC = data[:,1]   # [deg]
      self.Z = data[:,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:,3]   # [Mpc/h]
      self.coordY = data[:,4]   # [Mpc/h]
      self.coordZ = data[:,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:,6]   # [Mpc/h]
      self.dY = data[:,7]   # [Mpc/h]
      self.dZ = data[:,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:,10]   # [Mpc/h]
      self.dZKaiser = data[:,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:,12]   #[km/s]
      self.vY = data[:,13]   #[km/s]
      self.vZ = data[:,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:,15]  # [km/s]   from spherical catalo
      self.vTheta = data[:,16]   # [km/s]
      self.vPhi = data[:,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:,18]   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = data[:,19]
      self.Mvir = data[:,20]  # [M_sun]
      #
      # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      self.integratedTau = data[:,21]   # [dimless]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T v/c Tcmb
      self.integratedKSZ = data[:, 22] # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      self.integratedY = data[:, 23] # [sr]


   ##################################################################################
   ##################################################################################
   
   def addCatalog(self, newCat, save=False):
      """Combines the current catalog with a new catalog newCat.
      """
      # number of objects
      self.nObj += newCat.nObj
      #
      # sky coordinates and redshift
      self.RA = np.concatenate((self.RA, newCat.RA)) # [deg]
      self.DEC = np.concatenate((self.DEC, newCat.DEC))   # [deg]
      self.Z = np.concatenate((self.Z, newCat.Z))
      #
      # observed cartesian coordinates
      self.coordX = np.concatenate((self.coordX, newCat.coordX))   # [Mpc/h]
      self.coordY = np.concatenate((self.coordY, newCat.coordY))   # [Mpc/h]
      self.coordZ = np.concatenate((self.coordZ, newCat.coordZ))   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = np.concatenate((self.dX, newCat.dX))  # [Mpc/h]
      self.dY = np.concatenate((self.dY, newCat.dY))   # [Mpc/h]
      self.dZ = np.concatenate((self.dZ, newCat.dZ))  # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = np.concatenate((self.dXKaiser, newCat.dXKaiser))  # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = np.concatenate((self.dYKaiser, newCat.dYKaiser))   # [Mpc/h]
      self.dZKaiser = np.concatenate((self.dZKaiser, newCat.dZKaiser))   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = np.concatenate((self.vX, newCat.vX))   #[km/s]
      self.vY = np.concatenate((self.vY, newCat.vY))   #[km/s]
      self.vZ = np.concatenate((self.vZ, newCat.vZ))  #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = np.concatenate((self.vR, newCat.vR))  # [km/s]   from spherical catalo
      self.vTheta = np.concatenate((self.vTheta, newCat.vTheta))   # [km/s]
      self.vPhi = np.concatenate((self.vPhi, newCat.vPhi))  # [km/s]
      #
      # Stellar masses
      self.Mstellar = np.concatenate((self.Mstellar, newCat.Mstellar))   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = np.concatenate((self.hasM, newCat.hasM))
      self.Mvir = np.concatenate((self.Mvir, newCat.Mvir))  # [M_sun]
      #
      # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      self.integratedTau = np.concatenate((self.integratedTau, newCat.integratedTau))   # [dimless]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T v/c Tcmb
      self.integratedKSZ = np.concatenate((self.integratedKSZ, newCat.integratedKSZ)) # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      self.integratedY = np.concatenate((self.integratedY, newCat.integratedY)) # [sr]

      # Write the full catalog to the output path, if needed
      if save:
         self.writeCatalog()



   ##################################################################################
   ##################################################################################
   
   def plotFootprint(self):
      """Overlay a scatter plot of the catalog positions on top of a healpix map,
      here the AdvACT hit count map.
      """
      fig=plt.figure(0)
      #
      # hit count map for AdvACT
      path = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/healpix_f150_daynight_all_div_mono.fits"
      hHitMap = hp.read_map(path)
      hp.mollview(np.log(np.abs(hHitMap)+1.e-5), fig=0, title="", coord=None, cbar=False, unit='')
      #
      # scatter plot of the catalog
      hp.projscatter(self.RA, self.DEC, alpha=0.01, lonlat=True, marker='.', c='r', rasterized=True)
      #
      fig.savefig(self.pathFig+"/footprint_"+self.name+".pdf", dpi=1200)
      fig.clf()


   ##################################################################################
   ##################################################################################
   
   def printProperties(self):
      print "Catalog: "+self.nameLong
      print "Number of objects = "+str(self.nObj)
      print "with mass: "+str(np.sum(self.hasM))+", i.e. fraction "+str(np.sum(self.hasM)/self.nObj)
      print "Z: mean = "+str(np.mean(self.Z))+", median = "+str(np.median(self.Z))
      m = self.Mstellar[self.hasM==1]
      print "M_star [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
      m = self.Mvir[self.hasM==1]
      print "M_vir [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
   


   ##################################################################################
   ##################################################################################

   def compareV1dRms(self):
      """expected std dev of velocities for LCDM:
      comparison between using median z or z-distribution
      """
      # interpolate RMS 1d velocity for speed
      f = lambda z: self.U.v3dRms(0., z, W3d_sth) / np.sqrt(3.)
      Z = np.linspace(0., 1., 201)
      V1dRms = np.array(map(f, Z))
      f = interp1d(Z, V1dRms, kind='linear', bounds_error=False, fill_value='extrapolate')
      
      print "Expected v1d_rms = "+str(np.mean(np.array(map(f, self.Z))))+" km/s"
      print "Expected v1d_rms(z=z_mean) = "+str(f(np.mean(self.Z)))+" km/s"
      print "RMS v_r, v_theta, v_phi = "+str(np.std(self.vR))+", "+str(np.std(self.vTheta))+", "+str(np.std(self.vPhi))+" km/s"



   ##################################################################################
   ##################################################################################


   def plotHistograms(self):
      z0 = np.mean(self.Z)
      s2v1d = self.U.v3dRms(0., z0, W3d_sth)**2 / 3.
      
      # redshifts
      path = self.pathFig+"/hist_z.pdf"
      myHistogram(self.Z, nBins=71, lim=(0., 1.), path=path, nameLatex=r'$z$', semilogy=True)
      
      # spherical velocities
      path = self.pathFig+"/hist_vr.pdf"
      myHistogram(self.vR, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_r$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vtheta.pdf"
      myHistogram(self.vTheta, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_\theta$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vphi.pdf"
      myHistogram(self.vPhi, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_\phi$ [km/s]', doGauss=True)
      
      # cartesian velocities
      path = self.pathFig+"/hist_vx.pdf"
      myHistogram(self.vX, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_x$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vy.pdf"
      myHistogram(self.vY, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_y$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vz.pdf"
      myHistogram(self.vZ, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_z$ [km/s]', doGauss=True)
      
      # stellar masses
      path = self.pathFig+"/hist_mstellar.pdf"
      myHistogram(self.Mstellar, nBins=71, path=path, nameLatex=r'$M_\star$ [M$_\odot$]', semilogx=True, semilogy=True)

      # virial masses
      path = self.pathFig+"/hist_mvir.pdf"
      myHistogram(self.Mvir, nBins=71, path=path, nameLatex=r'$M_\text{vir}$ [M$_\odot$]', semilogx=True, semilogy=True)
      
      # comoving virial radius
      # need masses in Msun/h
      Par = zip(self.Mvir*self.U.bg.h, self.Z)
      f = lambda par: self.U.frvir(par[0], par[1])   # in: Msun/h, out: Mpc/h
      Rvir = np.array(map(f, Par))  # in Mpc/h
      #Rvir /= self.U.bg.h  # Mpc
      path = self.pathFig+"/hist_rvir.pdf"
      myHistogram(Rvir/self.U.bg.h, nBins=71, path=path, nameLatex=r'$R_\text{vir}$ [Mpc]', semilogx=True, semilogy=True)
      
      # virial angular radius
      Chi = np.array(map(self.U.bg.comoving_distance, self.Z)) # [Mpc/h]
      Thetavir = Rvir / Chi   # [rad]
      path = self.pathFig+"/hist_thetavir.pdf"
      x = Thetavir * (180.*60./np.pi)  # [arcmin]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\theta_\text{vir}$ [arcmin]', semilogx=True, semilogy=True)
      
      # integrated tau [arcmin^2]
      path = self.pathFig+"/hist_integratedtau.pdf"
      x = self.integratedTau * (180.*60./np.pi)**2 # [arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \tau$ [arcmin$^2$]', semilogx=True, semilogy=True)
      
      # mean tau within Rvir [dimless]
      path = self.pathFig+"/hist_meantauvir.pdf"
      x = self.integratedTau / (np.pi * Thetavir**2) # [dimless]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \tau / \left( \pi \theta_\text{vir} \right)$ [dimless]', semilogx=True, semilogy=True)


      # expected kSZ [muK*arcmin^2]
      path = self.pathFig+"/hist_ksz.pdf"
      x = self.integratedKSZ * (180.*60./np.pi)**2 # [muK*arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \delta T_\text{kSZ}$ [$\mu$K.arcmin$^2$]', doGauss=True, semilogy=True)

      # mean kSZ within Rvir [muK]
      path = self.pathFig+"/hist_meankszvir.pdf"
      x = self.integratedKSZ / (np.pi * Thetavir**2) # [muK]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \delta T_\text{kSZ} / \left( \pi \theta_\text{vir} \right)$ [$\mu$K]', doGauss=True, semilogy=True)

      # expected Y [arcmin^2]
      path = self.pathFig+"/hist_y.pdf"
      x = self.integratedY * (180.*60./np.pi)**2 # [arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; y_\text{tSZ}$ [arcmin$^2$]', semilogx=True, semilogy=True)

      # mean Y within Rvir [dimless]
      path = self.pathFig+"/hist_meanyvir.pdf"
      x = self.integratedY / (np.pi * Thetavir**2) # [dimless]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; y_\text{tSZ} / \left( \pi \theta_\text{vir} \right)$ [dimless]', semilogx=True, semilogy=True)

      # displacements?


##################################################################################
##################################################################################
