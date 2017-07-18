C**** Z1QX1N.BS1.F    Create Z1QX1N, Nonfractional ocean    2011/11/15
C****
C**** Fortran Compile and Go:  FCG Z1QX1N.BS1.CRE HNTR4.o
C****
C**** Z1QX1N.BS1.CRE creates a DataFile of observed surface fractions,
C**** topography, and water and ice thicknesses.
C****
C**** Check the following:
C****   "1QX1"
C****   "IM  ="
C****   "JM  ="
C****   "dLATM =" assumes equal spacing in latitude including poles
C****   "HNTR40"
C****   "FILOUT"
C****   Lines following "Set following grid cells to be continent"
C****   Lines following "Set following grid cells to be ocean"
C****
C**** Input files:
C****  Z2MX2M.NGDC = FOCEN2: Ocean Fraction (0 or 1)
C****                ZETOP2: Solid Topography (m) except for ice shelves
C****     Z10MX10M = FLAKES: Lake Fraction (0:1)
C****      ZICEHXH = dZGICH: Glacial Ice Thickness (m)
C****                FGICEH: Glacial Ice Fraction (0:1)
C****                ZSOLDH: Ice Top Topography (m)
C****       ZNGDC1 = FCONT1: Continent Fraction (0:1)
C****                FGICE1: Glacial Ice Fraction (0:1)
C****
      Integer*4,Parameter ::
     *  IM2 = 10800, JM2 = 5400,  !  2x2 minutes
     *  IMS =  2160, JMS = 1080,  !  10x10 minutes = 1/6 x 1/6 degrees
     *  IMH =   720, JMH =  360,  !  1/2 x 1/2 degrees
     *  IM1 =   360, JM1 =  180,  !  1x1 degrees
     *  IM  =   288, JM  =  180   !  1Qx1 degrees
      Real*4,Parameter :: dLATM = 180*60/JM  !  latitude spacing (minute
      Real*8,Parameter :: TWOPI = 6.283185307178586477d0,
     *                    AREAG = 2*TWOPI
      Real*4       FONLAT(JM),FOFLAT(JM)
      Real*8       SUMDFR,SUMDF,CONSTK
      Character    CHEM(IM)
      Character*80 IC,FILEIN,TITLEI, TITLE(13),TITLEF,FILOUT,FILEOF
      Common /Z2/ FOCEN2(IM2,JM2),ZETOP2(IM2,JM2),
     *            FCONT2(IM2,JM2),FGICE2(IM2,JM2),dZGIC2(IM2,JM2),
     *            ZSOLD2(IM2,JM2),ZSOLG2(IM2,JM2)
      Common /ZS/ FLAKES(IMS,JMS)
      Common /ZH/ dZGICH(IMH,JMH),FGICEH(IMH,JMH),ZSOLDH(IMH,JMH),
     *            FCON1H(IMH,JMH),FGIC1H(IMH,JMH),RGIC1H(IMH,JMH)
      Common /Z1/ FCONT1(IM1,JM1),FGICE1(IM1,JM1)
      Common /ZC/ FOCENF(IM,JM),FOCEAN(IM,JM),FLAKE(IM,JM),FLAND(IM,JM),
     *            dZOCEN(IM,JM),dZLAKE(IM,JM),FGRND(IM,JM),FGICE(IM,JM),
     *            dZGICE(IM,JM),ZSOLDG(IM,JM),ZATMO(IM,JM),
     *             ZSGLO(IM,JM), ZLAKE(IM,JM),ZGRND(IM,JM),ZSGHI(IM,JM)
      Real*8          DXYP,DXYPH,DXYP2
      Common /GEOMCB/ DXYP(JM), DXYPH(JMH)
      Common /WORKCB/ WT2(IM2,JM2)
C****
C**** Output records:
C****
      Data TITLE /
     1'FOCEAN: 0 or 1, Bering Strait 1 cell wide         GISS 1Qx1',
     2'FLAKE: Lake Surface Fraction (0:1)                GISS 1Qx1',
     3'FGRND: Ground Surface Fraction (0:1)              GISS 1Qx1',
     4'FGICE: Glacial Ice Surface Fraction (0:1)         GISS 1Qx1',
     5'ZATMO: Atmospheric Topography (m)                 ETOPO2 1Qx1',
     6'dZOCEN: Ocean Thickness (m)                       ETOPO2 1Qx1',
     7'dZLAKE: Lake Thickness (m)                        ETOPO2 1Qx1',
     8'dZGICE: Glacial Ice Thickness (m)                 Ekholm,Bamber',
     9'ZSOLDG: Solid Ground Topography (m)               ETOPO2 1Qx1',
     A'ZSGLO: Lowest Solid Topography (m)                ETOPO2 1Qx1',
     B'ZLAKE: Lake Surface Topography (m)                ETOPO2 1Qx1',
     C'ZGRND: Topography Break between Ground and GIce   ETOPO2 1Qx1',
     D'ZSGHI: Highest Solid Topography (m)               ETOPO2 1Qx1'/
      Data TITLEF /
     1'FOCENF: Fractional Ocean Cover (0:1)              ETOPO2 1Qx1'/
C****
      FILOUT = 'Z1QX1N.BS1'
      FILEOF = 'Z1QX1F.BS1'  !  Fractional ocean cover
      Call SGEOM
C****
C**** Read input DataFiles
C****
      Call GetEnv ('IC',IC)
C**** Read in Z2MX2M.NGDC
      FILEIN = Trim(IC) // '/Z/Z2MX2M.NGDC'
      Open  (1, File=FILEIN, Form='Unformatted', Status='Old',
     &    Convert='BIG_ENDIAN')
      Read  (1) TITLEI,FOCEN2
      Write (6,*) 'FOCEN2 read from Z2MX2M.NGDC: ',TITLEI
      Read  (1) TITLEI,ZETOP2
      Write (6,*) 'ZETOP2 read from Z2MX2M.NGDC: ',TITLEI
      Close (1)
C**** Read in Z10MX10M
      FILEIN = Trim(IC) // '/Z/Z10MX10M'
      Open  (2, File=FILEIN, Form='Unformatted', Status='Old',
     &    Convert='BIG_ENDIAN')
      Read  (2)
      Read  (2) TITLEI,FLAKES
      Write (6,*) 'FLAKES read from Z10MX10M: ',TITLEI
      Close (2)
C**** Read in ZICEHXH
      FILEIN = Trim(IC) // '/Z/ZICEHXH'
      Open  (3, File=FILEIN, Form='Unformatted', Status='Old',
     &    Convert='BIG_ENDIAN')
      Read  (3) TITLEI,dZGICH
      Write (6,*) 'dZGICH read from ZICEHXH: ',TITLEI
      Read  (3) TITLEI,FGICEH
      Write (6,*) 'FGICEH read from ZICEHXH: ',TITLEI
      Read  (3) TITLEI,ZSOLDH
      Write (6,*) 'ZSOLDH read from ZICEHXH: ',TITLEI
      Close (3)
C**** Read in ZNGDC1
      FILEIN = Trim(IC) // '/Z/ZNGDC1'
      Open  (4, File=FILEIN, Form='Unformatted', Status='Old',
     &    Convert='BIG_ENDIAN')
      Read  (4)
      Read  (4)
      Read  (4)
      Read  (4) TITLEI,FCONT1
      Write (6,*) 'FCONT1 read from ZNDGC1: ',TITLEI
      Read  (4)
      Read  (4) TITLEI,FGICE1
      Write (6,*) 'FGICE1 read from ZNGDC1: ',TITLEI
      Close (4)
      Write (0,*) 'Input files read in.'
C****
C$OMP ParallelDo   Private (J2)
      Do 10 J2=1,JM2
   10 WT2(:,J2) = 1
C****
C**** Add small ice cap and glacier data to FGICEH and dZGICH
C**** north of Antarctic area.
C**** Continental cells north of 78N are entirely glacial ice.
C****
      FGICE1(:,JM1*14/15+1:JM1) = FCONT1(:,JM1*14/15+1:JM1)
      Call HNTR40 (IM1,JM1,0.,60., IMH,JMH,0.,30., 0.)
      Call HNTR4  (WT2,FCONT1, FCON1H)
      Call HNTR4  (WT2,FGICE1, FGIC1H)
C**** RGIC1H = areal ratio of glacial ice to continent
C**** For smaller ice caps and glaciers, dZGICH = CONSTK * RGIC1H^.3
C**** Constant is chosen so that average value of dZGICH is 264.7 m
C**** 264.7  =  sum(DXYP*FGIC1H*dZGICH) / sum(DXYP*FGIC1H)  =
C****        =  CONSTK * sum(DXYP*FGIC1H*RGIC1H^.3) / sum(DXYP*FGIC1H)
      SUMDFR = 0  ;  SUMDF = 0
      Do 30 JH=1+JMH/6,JMH
      Do 20 IH=1,IMH
      If (FGICEH(IH,JH) > 0)  FGIC1H(IH,JH) = 0  !  ignore Greenland
   20 RGIC1H(IH,JH) = FGIC1H(IH,JH) / (FCON1H(IH,JH)+1e-20)
      SUMDFR = SUMDFR + DXYPH(JH)*Sum(FGIC1H(:,JH)*RGIC1H(:,JH)**.3)
   30 SUMDF  = SUMDF  + DXYPH(JH)*Sum(FGIC1H(:,JH))
      CONSTK = 264.7 * SUMDF / SUMDFR
C**** Replace FGICEH and dZGICH away from Greenland
      Do 40 JH=JMH/6+1,JMH
      Do 40 IH=1,IMH
      If (FGICEH(IH,JH) == 0)  Then
          FGICEH(IH,JH) = FGIC1H(IH,JH)
          dZGICH(IH,JH) = CONSTK * RGIC1H(IH,JH)**.3  ;  EndIf
   40 Continue
C****
C**** ETOPO2 treats Antarctic ice shelves as ocean.
C**** When this happens ETOPO2 data is replaced with interpolated data
C**** from FGICEH and dZGICH.  Resulting files are:
C**** FOCEN2 = Ocean fraction (0 or 1) correct for Antarctic ice shelves
C**** FCONT2 = Continental fraction (0 or 1)
C**** FGICE2 = Glacial ice fraction (0 or 1)
C**** dZGIC2 = Thickness of glacial ice (m)
C**** ZSOLD2 = Solid topography (m)        (above ice)
C**** ZSOLG2 = Solid ground topography (m) (beneath ice)
C****
      Call HNTR40 (IMH,JMH,0.,30., IM2,JM2,0.,2., 0.)
      Call HNTR4  (WT2   ,FGICEH, FGICE2)
      Call HNTR4  (FGICEH,dZGICH, dZGIC2)
      Call HNTR4  (FGICEH,ZSOLDH, ZSOLD2)
C**** North of Antarctic area: 60S to 90N
C$OMP ParallelDo   Private (J2)
      Do 60 J2=JM2/6+1,JM2
      FCONT2(:,J2) = 1 - FOCEN2(:,J2)
      FGICE2(:,J2) = FGICE2(:,J2) * FCONT2(:,J2)
      dZGIC2(:,J2) = dZGIC2(:,J2) * FCONT2(:,J2)
      ZSOLD2(:,J2) = ZETOP2(:,J2)
   60 ZSOLG2(:,J2) = ZETOP2(:,J2) - dZGIC2(:,J2)
C**** Antarctic area: 90S to 60S
C$OMP ParallelDo   Private (I2,J2)
      Do 70 J2=1,JM2/6
      Do 70 I2=1,IM2
      If (FOCEN2(I2,J2) == 0)  Then
C         FOCEN2(I2,J2) = 0
          FCONT2(I2,J2) = 1
          FGICE2(I2,J2) = 1
C         dZGIC2(I2,J2) = dZGIC2(I2,J2)  !  ZETOP2 has 2m and other low
C         ZSOLD2(I2,J2) = ZSOLD2(I2,J2)       values over ice shelves
            If (ZETOP2(I2,J2) >= 100)  ZSOLD2(I2,J2) = ZETOP2(I2,J2)
          ZSOLG2(I2,J2) = ZSOLD2(I2,J2) - dZGIC2(I2,J2)
          GoTo 70  ;  EndIf
      If (FGICE2(I2,J2) <= .5)  Then  !  and FOCEN2(I2,J2) == 1
C         FOCEN2(I2,J2) = 1
          FCONT2(I2,J2) = 0
          FGICE2(I2,J2) = 0
          dZGIC2(I2,J2) = 0
          ZSOLD2(I2,J2) = ZETOP2(I2,J2)
          ZSOLG2(I2,J2) = ZETOP2(I2,J2)
          GoTo 70  ;  EndIf
      If (FGICE2(I2,J2) > .5)  Then  !  and FOCEN2(I2,J2) == 1
          FOCEN2(I2,J2) = 0
          FCONT2(I2,J2) = 1
          FGICE2(I2,J2) = 1
          dZGIC2(I2,J2) = ZSOLD2(I2,J2) - ZETOP2(I2,J2)
C         ZSOLD2(I2,J2) = ZSOLD2(I2,J2)
          ZSOLG2(I2,J2) = ZETOP2(I2,J2)  ;  EndIf
   70 Continue
C****
C**** FOCEAN: Ocean Surface Fraction (0:1)
C****
C**** Fractional ocean cover FOCENF is interpolated from FOAAH2
      Call HNTR40 (IM2,JM2,0.,2., IM,JM,0.,dLATM, 0.)
      Call HNTR4P (WT2,FOCEN2, FOCENF)
C**** FOCEAN (0 or 1) is rounded from FOCEAN
      FOCEAN(:,:) = Nint (FOCENF(:,:))
C**** Set following grid cells to be continent
      FOCEAN( 84, 18) = 0
      FOCEAN( 85, 18) = 0
      FOCEAN(236, 82) = 0
      FOCEAN(242, 82) = 0
      FOCEAN(245, 82) = 0
      FOCEAN(224,101) = 0
      FOCEAN( 53,119) = 0
      FOCEAN(171,125) = 0  !  Cyprus
      FOCEAN(164,126) = 0  !  Crete
      FOCEAN(158,129) = 0
      FOCEAN(158,130) = 0
      FOCEAN(242,131) = 0
      FOCEAN(263,136) = 0
      FOCEAN(258,137) = 0
      FOCEAN(258,138) = 0
      FOCEAN( 46,139) = 0
      FOCEAN(258,139) = 0
      FOCEAN(275,152) = 0
      FOCEAN(  8,156) = 0
      FOCEAN( 10,156) = 0
      FOCEAN( 12,157) = 0
      FOCEAN(172,157) = 0
      FOCEAN(202,157) = 0
      FOCEAN( 69,159) = 0
      FOCEAN(204,159) = 0
      FOCEAN( 62,167) = 0
      FOCEAN( 73,171) = 0
      FOCEAN( 75,171) = 0
      FOCEAN( 78,171) = 0
C**** Set following grid cells to be ocean
      FOCEAN(179,105) = 1
      FOCEAN( 54,119) = 1
      FOCEAN(241,131) = 1
      FOCEAN(258,143) = 1
      FOCEAN(165,150) = 1
      FOCEAN(274,152) = 1
      FOCEAN( 15,154) = 1
      FOCEAN( 92,155) = 1
      FOCEAN( 13,157) = 1
      FOCEAN(173,157) = 1
      FOCEAN(176,157) = 1
      FOCEAN(203,157) = 1
      FOCEAN( 55,159) = 1
      FOCEAN(103,159) = 1
      FOCEAN(203,159) = 1
      FOCEAN( 67,160) = 1
      FOCEAN( 68,160) = 1
      FOCEAN( 79,160) = 1
      FOCEAN(199,160) = 1
      FOCEAN(126,161) = 1
      FOCEAN( 68,162) = 1
      FOCEAN( 75,165) = 1
      FOCEAN(225,169) = 1
C**** Average non-fractional and fractional ocean covers over latitude
      Write (6,910)
      Do 110 J=JM,1,-1
      FOFLAT(J) = Sum(FOCENF(:,J))
      FONLAT(J) = Sum(FOCEAN(:,J))
  110 Write (6,911) J,FOFLAT(J),FONLAT(J),FOFLAT(J)-FONLAT(J)
      FOFSH = IM*JM * Sum(FOFLAT(1:JM/2)*DXYP(1:JM/2)) / AREAG
      FONSH = IM*JM * Sum(FONLAT(1:JM/2)*DXYP(1:JM/2)) / AREAG
      FOFNH = IM*JM * Sum(FOFLAT(1+JM/2:JM)*DXYP(1+JM/2:JM)) / AREAG
      FONNH = IM*JM * Sum(FONLAT(1+JM/2:JM)*DXYP(1+JM/2:JM)) / AREAG
      Write (6,912) 'NH',FOFNH,FONNH,FOFNH-FONNH
      Write (6,912) 'SH',FOFSH,FONSH,FOFSH-FONSH
C**** Write fractional ocean cover to output file FILEOF
      Open  (10, File=FILEOF, Form='Unformatted')
      Write (10) TITLEF,FOCENF
      Close (10)
      Write (6,*) 'FOCENF written to unit 10: ',TITLEF
C****
C**** FLAKE: Lake Surface Fraction (0:1)
C****
C**** FLAKE is interpolated from FLAKES
      Call HNTR40 (IMS,JMS,0.,10., IM,JM,0.,dLATM, 0.)
      Call HNTR4P (WT2,FLAKES, FLAKE)
C**** Antarctica and Arctic area have no lakes
      FLAKE(:,1:JM/6) = 0             !  90:60 S
      FLAKE(:,JM*14/15+1:JM) = 0      !  78:90 N
      FLAKE(:IM/2,JM*41/45+1:JM) = 0  !  74:90 N, 0:180 W
      Do 160 J=JM*5/6,JM*11/12                   !  southern
      Do 160 I=IM/3+1,Nint(.75*IM*(J-JM*.3)/JM)  !  Greenland
  160 FLAKE(I,J) = 0
C**** Apportion FLAKE to the nonocean fraction and round to 1/256
      FLAKE(:,:) = FLAKE(:,:)*(1-FOCEAN(:,:)) / (1-FOCENF(:,:)+1e-20)
      FLAKE(:,:) = Nint(FLAKE(:,:)*256) / 256.
C****
C**** FGICE: Glacial Ice Surface Fraction (0:1)
C****
C**** FGICE is interpolated from FGICE2
      Call HNTR40 (IM2,JM2,0.,2., IM,JM,0.,dLATM, 0.)
      Call HNTR4P (FCONT2,FGICE2, FGICE)

      Open  (11, File='int_f', Form='Unformatted')
      Write (11) TITLE(1),FOCEAN
      write (11) TITLE(2),FLAKE
      write (11) TITLE(4),FGICE
      close (11)


C**** Antarctica is entirely glacial ice, no lakes nor ground
      FGICE(:,1:JM/6) = 1 - FOCEAN(:,1:JM/6)
C**** Continental cells north of 78N are entirely glacial ice
      FGICE(:,JM*14/15+1:JM) = 1 - FOCEAN(:,JM*14/15+1:JM)
C**** There is no glacial ice over oceans
      FGICE(:,JM/6+1:JM) = FGICE(:,JM/6+1:JM) * (1-FOCEAN(:,JM/6+1:JM))
C**** Round FGICE to nearest 1/256
      FGICE(:,:) = Nint(FGICE(:,:)*256) / 256.


C**** Check that FGICE is between 0 and 1
C**** If FGICE+FLAKE exceeds 1, reduce FLAKE
      Do 210 J=JM/6+1,JM
      Do 210 I=1,IM
      If (FGICE(I,J) < 0)  Then
          Write (0,*) '210: FGICE < 0:',I,J,FGICE(I,J)
          FGICE(I,J) = 0  ;  EndIf
      If (FGICE(I,J) > 1)  Then
          Write (0,*) '210: FGICE > 1:',I,J,FGICE(I,J)
          FGICE(I,J) = 1  ;  EndIf
      If (FLAKE(I,J)+FGICE(I,J)+FOCEAN(I,J) > 1)  Then
          Write (0,*) '210: FGICE+FLAKE+FOCEAN > 1:',I,J,FGICE(I,J),
     *          FLAKE(I,J),FOCEAN(I,J)
          FLAKE(I,J) = 1 - FGICE(I,J) - FOCEAN(I,J)  ;  EndIf
  210 Continue
C**** Replace land cells without vegetation with glacial ice in Z1QX1N
c     FGICE(35,52) = 1 - FLAKE(35,52)     4x3 Model
C****
C**** FGRND: Surface Fraction of Ground (0:1)
C****
      FGRND(:,:) = 1 - FOCEAN(:,:) - FLAKE(:,:) - FGICE(:,:)
C**** Check that FGRND is between 0 and 1
      Do 260 J=1,JM
      Do 260 I=1,IM
  260 If (FGRND(I,J) < 0 .or. FGRND(I,J) > 1)  Write (0,*)
     *  'Error: FGRND=',I,J,FGRND(I,J),FOCEAN(I,J),FLAKE(I,J),FGICE(I,J)
C****
C**** dZOCEN: Ocean Thickness (m)
C****
C     Call HNTR40 (IM2,JM2,0.,2., IM,JM,0.,dLATM, 0.)
      Call HNTR4P (FOCEN2,ZSOLG2, dZOCEN)
      dZOCEN(:,:) = - dZOCEN(:,:)*FOCEAN(:,:)
C**** Check that dZOCEN is positive
      Do 310 J=1,JM
      Do 310 I=1,IM
      If (FOCEAN(I,J) == 1 .and. dZOCEN(I,J) <= 0)
     *  Write (0,*) 'Error: dZOCEN <= 0 at:',I,J,dZOCEN(I,J)
  310 Continue
C****
C**** dZGICE: Glacial Ice Thickness (m)
C****
C     Call HNTR40 (IM2,JM2,0.,2., IM,JM,0.,dLATM, 0.)
      Call HNTR4P (FGICE2,dZGIC2, dZGICE)
      Do 360 J=1,JM
      Do 360 I=1,IM
      If (FGICE(I,J) > 0)
     *   Then  ;  dZGICE(I,J) = Max (dZGICE(I,J), 1d0)
         Else  ;  dZGICE(I,J) = 0  ;  EndIf
  360 Continue
C****
C**** ZATMO  = Atmospheric topography (m)
C**** dZLAKE = Mean lake thickness (m)
C**** ZSOLDG = Solid ground topography (m)
C**** ZSGLO  = Lowest value of ZSOLD2 in model cell (m)
C**** ZLAKE  = Surface lake topography (m)
C**** ZGRND  = Altitude break between ground and land ice (m)
C**** ZSGHI  = Highest value of ZSOLD2 in model cell (m)
C****
      ZSOLDG(:,:) = - dZOCEN(:,:)  !  solid ground topography of ocean
      Call Z (FOCEN2,ZSOLD2,ZSOLG2, FOCEAN,FLAKE,FGRND,
     *        ZATMO,dZLAKE,ZSOLDG,ZSGLO,ZLAKE,ZGRND,ZSGHI)
C****
C**** Reset ZATMO, dZOCEN and ZLAKE by hand if FLAKE(I,J) == 1
C****
C**** Caspian Sea
      dZLAKE(186,128) = dZLAKE(186,128) - 30 - ZLAKE(186,128)
      dZLAKE(187,128) = dZLAKE(187,128) - 30 - ZLAKE(187,128)
      dZLAKE(185,129) = dZLAKE(185,129) - 30 - ZLAKE(185,129)
      dZLAKE(186,129) = dZLAKE(186,129) - 30 - ZLAKE(186,129)
      dZLAKE(185,130) = dZLAKE(185,130) - 30 - ZLAKE(185,130)
      dZLAKE(186,130) = dZLAKE(186,130) - 30 - ZLAKE(186,130)
      dZLAKE(186,131) = dZLAKE(186,131) - 30 - ZLAKE(186,131)
      dZLAKE(185,132) = dZLAKE(185,132) - 30 - ZLAKE(185,132)
      dZLAKE(184,133) = dZLAKE(184,133) - 30 - ZLAKE(184,133)
      dZLAKE(185,133) = dZLAKE(185,133) - 30 - ZLAKE(185,133)
      dZLAKE(184,134) = dZLAKE(184,134) - 30 - ZLAKE(184,134)
      dZLAKE(183,135) = dZLAKE(183,135) - 30 - ZLAKE(183,135)
      dZLAKE(184,135) = dZLAKE(184,135) - 30 - ZLAKE(184,135)
       ZATMO(186,128) = - 30  ;  ZLAKE(186,128) = -30
       ZATMO(187,128) = - 30  ;  ZLAKE(187,128) = -30
       ZATMO(185,129) = - 30  ;  ZLAKE(185,129) = -30
       ZATMO(186,129) = - 30  ;  ZLAKE(186,129) = -30
       ZATMO(185,130) = - 30  ;  ZLAKE(185,130) = -30
       ZATMO(186,130) = - 30  ;  ZLAKE(186,130) = -30
       ZATMO(186,131) = - 30  ;  ZLAKE(186,131) = -30
       ZATMO(185,132) = - 30  ;  ZLAKE(185,132) = -30
       ZATMO(184,133) = - 30  ;  ZLAKE(184,133) = -30
       ZATMO(185,133) = - 30  ;  ZLAKE(185,133) = -30
       ZATMO(184,134) = - 30  ;  ZLAKE(184,134) = -30
       ZATMO(183,135) = - 30  ;  ZLAKE(183,135) = -30
       ZATMO(184,135) = - 30  ;  ZLAKE(184,135) = -30
C**** Aral Sea
      dZLAKE(192,135) = dZLAKE(192,135) + 53 - ZLAKE(192,135)
       ZATMO(192,135) = 53  ;  ZLAKE(192,135) = 53
C**** Lake Superior
      dZLAKE(75,138) = dZLAKE(75,138) + 602 - ZLAKE(75,138)
       ZATMO(75,138) = 602  ;  ZLAKE(75,138) = 602
C****
      Do 410 J=1,JM
      Do 410 I=1,IM
  410 If (FLAKE(I,J) == 1)  Write (0,*) 'FLAKE == 1 at I,J=',I,J,
     *                        ZATMO(I,J),dZLAKE(I,J),ZLAKE(I,J)
C****
C**** Write output files to disk
C****
      Open  (11, File=FILOUT, Form='Unformatted')
      Write (11) TITLE(1),FOCEAN
      Write (6,*) 'FOCEAN written on unit 11: ',TITLE(1)
      Write (11) TITLE(2),FLAKE
      Write (6,*) 'FLAKE written on unit 11: ',TITLE(2)
      Write (11) TITLE(3),FGRND
      Write (6,*) 'FGRND written on unit 11: ',TITLE(3)
      Write (11) TITLE(4),FGICE
      Write (6,*) 'FGICE written on unit 11: ',TITLE(4)
      Write (11) TITLE(5),ZATMO
      Write (6,*) 'ZATMO written on unit 11: ',TITLE(5)
      Write (11) TITLE(6),dZOCEN
      Write (6,*) 'dZOCEN written on unit 11: ',TITLE(6)
      Write (11) TITLE(7),dZLAKE
      Write (6,*) 'dZLAKE written on unit 11: ',TITLE(7)
      Write (11) TITLE(8),dZGICE
      Write (6,*) 'dZGICE written on unit 11: ',TITLE(8)
      Write (11) TITLE(9),ZSOLDG
      Write (6,*) 'ZSOLDG written on unit 11: ',TITLE(9)
      Write (11) TITLE(10),ZSGLO
      Write (6,*) 'ZSGLO written on unit 11: ',TITLE(10)
      Write (11) TITLE(11),ZLAKE
      Write (6,*) 'ZLAKE written on unit 11: ',TITLE(11)
      Write (11) TITLE(12),ZGRND
      Write (6,*) 'ZGRND written on unit 11: ',TITLE(12)
      Write (11) TITLE(13),ZSGHI
      Write (6,*) 'ZSGHI written on unit 11: ',TITLE(13)
      Close (11)
C****
  910 Format (/
     *' Comparison between Fractional and Non-fractional Ocean Cover' //
     *'         # of      # of     differ' /
     *'         fract    NOfrac      in #' /
     *'   J     cells     cells     cells' /
     *'   =     =====     =====     =====')
  911 Format (I4,3F10.2)
  912 Format (2X,A2,3F10.2)
      End

      Subroutine SGEOM
C****
C**** SGEOM calculates the sperical area of grid cells
C****
      Implicit  Real*8 (A-H,O-Z)
      Integer*4,Parameter ::
     *  IM2 = 10800, JM2 = 5400,  !  2x2 minutes
     *  IMH =   720, JMH =  360,  !  1/2 x 1/2 degrees
     *  IM  =   288, JM  =  180   !  1Qx1 degrees
      Real*8,Parameter :: TWOPI = 6.283185307178586477d0
      Common /GEOMCB/ DXYP(JM),DXYPH(JMH),DXYP2(JM2)
C****
C**** Area for IM,JM resolution
C****
      dLON =    TWOPI / IM
      dLAT = .5*TWOPI / JM
      Do 10 J=1,JM
      SINS = Sin (dLAT*(J-JM/2-1))
      SINN = Sin (dLAT*(J-JM/2))
   10 DXYP(J) = dLON * (SINN - SINS)
C****
C**** Area for IMH,JMH resolution
C****
      dLON =    TWOPI / IMH
      dLAT = .5*TWOPI / JMH
      Do 20 JH=1,JMH
      SINS = Sin (dLAT*(JH-JMH/2-1))
      SINN = Sin (dLAT*(JH-JMH/2))
   20 DXYPH(JH) = dLON * (SINN - SINS)
C****
C**** Area for IMH,JMH resolution
C****
      dLON =    TWOPI / IM2
      dLAT = .5*TWOPI / JM2
      Do 30 J2=1,JM2
      SINS = Sin (dLAT*(J2-JM2/2-1))
      SINN = Sin (dLAT*(J2-JM2/2))
   30 DXYP2(J2) = dLON * (SINN - SINS)
C****
      Return
      End

      Subroutine Z (FOCEN2,ZSOLD2,ZSOLG2, FOCEAN,FLAKE,FGRND,
     *              ZATMO,dZLAKE,ZSOLDG,ZSGLO,ZLAKE,ZGRND,ZSGHI)
C****
C**** Input:  FOCEN2 = ocean fraction at 2 x 2 (minute)
C****         ZSOLD2 = solid topography (above ice) at 2 x 2 (minute)
C****         ZSOLG2 = solid ground topography at 2 x 2 (minute)
C****
C**** Output: ZATMO  = atmospheric topography (m)
C****         dZLAKE = mean lake thickness (m)
C****         ZSOLDG = solid ground topography (m)
C****         ZSGLO  = lowest value of ZSOLD2 in model cell (m)
C****         ZLAKE  = surface lake topography (m)
C****         ZGRND  = altitude break between ground and land ice (m)
C****         ZSGHI  = highes value of ZSOLD2 in model cell (m)
C****
      Implicit  None
      Integer*4,Parameter ::
     *  IM2 = 10800, JM2 = 5400,  !  2x2 minutes
     *  IMH =   720, JMH =  360,  !  2x2 minutes
     *  IM  =   288, JM  =  180   !  1Qx1 degrees
      Real*4 FOCEN2(IM2,JM2),ZSOLD2(IM2,JM2),ZSOLG2(IM2,JM2)
      Real*4 FOCEAN(IM,JM),FLAKE(IM,JM),FGRND(IM,JM),ZATMO(IM,JM),
     *       dZLAKE(IM,JM),ZLAKE(IM,JM),ZGRND(IM,JM),
     *       ZSOLDG(IM,JM),ZSGLO(IM,JM),ZSGHI(IM,JM)
      Real*8 DXYP,DXYPH,DXYP2
      Common /GEOMCB/ DXYP(JM),DXYPH(JMH),DXYP2(JM2)
C**** Local variables
      Integer*4 I,J,IMAX, I2,I21,I2M,J2,J21,J2M, N,NM,NNEW, NLAKE,NGRND
      Real*4    ZN(1296000),AN(1296000)
      Real*8    SAREA,SAZSG, ANSUM,VNSUM, ZLBOT, AZATMO
C****
      If (IM2*JM2/JM > 1296000)  GoTo 800
C****
      Do 310 J=1,JM
      J21 = (J-1)*JM2/JM + 1  ;  J2M = J*JM2/JM
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 300 I=1,IMAX
      I21 = (I-1)*IM2/IM + 1  ;  I2M = I*IM2/IM  ;  If(IMAX==1) I2M=IM2

      If (FOCEAN(I,J) == 0)  GoTo 200
C****
C**** (I,J) is an ocean cell
C****
       ZATMO(I,J) = 0
      dZLAKE(I,J) = 0
C     ZSOLDG(I,J) = - dZOCEN(I,J)  !  already filled in
       ZLAKE(I,J) = 0
       ZGRND(I,J) = 0
       ZSGHI(I,J) = 0
       ZSGLO(I,J) = 999999
      Do 110 J2=J21,J2M
      Do 110 I2=I21,I2M
  110 If (ZSGLO(I,J) > ZSOLD2(I2,J2) .and. FOCEN2(I2,J2) == 1)
     *    ZSGLO(I,J) = ZSOLD2(I2,J2)
      If (ZSGLO(I,J) == 999999)  GoTo 811
      GoTo 300
C****
C**** (I,J) is a continental cell
C****
C**** Order 2-minute continental cells within (I,J) and sum their area
  200 SAREA = 0  ;  SAZSG = 0
      NM = 0



      Do 240 J2=J21,J2M
      Do 240 I2=I21,I2M

      If (FOCEN2(I2,J2) == 1)  GoTo 240
      Do 210 NNEW=1,NM
  210 If (ZSOLD2(I2,J2) < ZN(NNEW))  GoTo 220
      NNEW = NM+1
C**** Insert NNEW into the arrays 1:NM
  220 NM = NM+1
      Do 230 N=NM,NNEW+1,-1
      AN(N) = AN(N-1)
  230 ZN(N) = ZN(N-1)
      AN(NNEW) =  DXYP2(J2)
      ZN(NNEW) = ZSOLD2(I2,J2)
      SAREA = SAREA + DXYP2(J2)
      SAZSG = SAZSG + DXYP2(J2)*ZSOLG2(I2,J2)
  240 Continue

      If (SAREA == 0)  GoTo 824
C**** Determine ZSOLDG
      ZSOLDG(I,J) = SAZSG / SAREA
C**** Determine ZSGLO
      ZSGLO(I,J) = ZN(1)
C**** Determine ZLAKE and dZLAKE
      ANSUM = 0  ;  VNSUM = 0  !  accumulated area and volume
      Do 250 NLAKE=1,NM
      ANSUM = ANSUM + AN(NLAKE)
      VNSUM = VNSUM + AN(NLAKE)*ZN(NLAKE)
  250 If (ANSUM > SAREA*FLAKE(I,J))  GoTo 260
      NLAKE = NM
  260 ZLAKE(I,J) = ZN(NLAKE)
      If (FLAKE(I,J) > 0)  Then
         ZLBOT = (VNSUM - (ANSUM - SAREA*FLAKE(I,J))*ZN(NLAKE)) /
     /           (SAREA*FLAKE(I,J))
         dZLAKE(I,J) = Max (ZN(NLAKE)-ZLBOT, 1d0)
      Else
         dZLAKE(I,J) = 0  ;  EndIf
C**** Determine ZATMO
      AZATMO = ANSUM*ZLAKE(I,J)
      Do 270 N=NLAKE+1,NM
  270 AZATMO = AZATMO + AN(N)*ZN(N)
      ZATMO(I,J) = AZATMO / SAREA
C**** Determine ZGRND
      ANSUM = ANSUM - AN(NLAKE)
      Do 280 NGRND=NLAKE,NM
      ANSUM = ANSUM + AN(NGRND)
  280 If (ANSUM > SAREA*(FLAKE(I,J)+FGRND(I,J)))  GoTo 290
      NGRND = NM
  290 ZGRND(I,J) = ZN(NGRND)
C**** Determine ZSGHI
      ZSGHI(I,J) = ZN(NM)
  300 Continue
C**** Replicate Z data to all longitudes at poles
      If (J==1 .or. J==JM)  Then
         ZATMO (2:IM,J) = ZATMO (1,J)
         dZLAKE(2:IM,J) = dZLAKE(1,J)
         ZSOLDG(2:IM,J) = ZSOLDG(1,J)
         ZSGLO (2:IM,J) = ZSGLO (1,J)
         ZLAKE (2:IM,J) = ZLAKE (1,J)
         ZGRND (2:IM,J) = ZGRND (1,J)
         ZSGHI (2:IM,J) = ZSGHI (1,J)  ;  EndIf
  310 Continue
      Return
C****
  800 Write (0,*) 'IM2*JM2/IM*JM > 1296000.  Increase 1296000.'
      Stop 800
  811 Write (0,*) 'Ocean cell I,J has no ocean area on 2-minute grid.'
      Write (0,*) 'I,J,I21,I2M,J21,J2M=',I,J,I21,I2M,J21,J2M
      Stop 814
  824 Write (0,*)
     *'Continental cell I,J has no continental area on 2-minute grid.'
      Write (0,*) 'I,J,I21,I2M,J21,J2M=',I,J,I21,I2M,J21,J2M
      Stop 824
      End
