module help_hntr
use iso_c_binding

implicit none

contains

subroutine write_hntr40(IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS) bind(C)
    Integer :: IMA,JMA,IMB,JMB
    Real*4  :: OFFIA,DLATA, OFFIB,DLATB, DATMIS
    ! ---------- Local vars

      Integer :: INA,JNA,INB,JNB, IMIN,IMAX,JMIN,JMAX, &
                 IA,IB,JA,JB, IBp1
      Real*4  :: DATMCB
      Real*8  :: SINA,SINB, FMIN,FMAX,GMIN,GMAX, &
                 DIA,DIB, RIA,RIB,RJA,RJB, FJEQA,FJEQB
    Common /HNTRCB/ SINA(0:5401),SINB(0:5401), &
                      FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401), &
                      IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401), &
                      DATMCB, INA,JNA,INB,JNB

    call hntr40(ima, jma, offia, dlata, imb, jmb, offib, dlatb, datmis)

    open (1, File='hntr4_common', Form='Unformatted', Status='REPLACE')
    write (1) datmcb, ina, jna, inb, jnb
    write (1) sina, sinb, fmin, fmax, gmin, gmax, imin, imax, jmin, jmax
    close (1)
end subroutine write_hntr40

subroutine call_hntr4(WTA, lwta, A, la, B, lb) bind(C)
    integer :: lwta, la, lb
    Real*4  :: WTA(lwta), A(la), B(lb)

    call HNTR4(WTA, A, B)

end subroutine call_hntr4


subroutine write_sgeom() bind(C)
    real*8 dxyp, dxyph, dxyp2
      Integer*4,Parameter ::    &
       IM2 = 10800, JM2 = 5400, &  !  2x2 minutes
       IMH =   720, JMH =  360, &  !  1/2 x 1/2 degrees
       IM  =   288, JM  =  180   !  1Qx1 degrees
    Common /GEOMCB/ DXYP(JM),DXYPH(JMH),DXYP2(JM2)

    call sgeom

    open (1, File='sgeom_common', Form='Unformatted', Status='REPLACE')
    write (1) dxyp,dxyph,dxyp2
    close (1)
end subroutine write_sgeom


end module help_hntr
