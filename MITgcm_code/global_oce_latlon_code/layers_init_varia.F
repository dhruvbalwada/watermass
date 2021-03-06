C $Header: /u/gcmpack/MITgcm/pkg/layers/layers_init_varia.F,v 1.8 2014/07/08 19:02:42 jmc Exp $
C $Name: checkpoint65c $

#include "LAYERS_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      SUBROUTINE LAYERS_INIT_VARIA( myThid )

C ===================================================================
C     Initialize LAYERS variables.
C ===================================================================

      IMPLICIT NONE
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "LAYERS_SIZE.h"
#include "LAYERS.h"

C  INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid

#ifdef ALLOW_LAYERS

C     === Local variables ===
      INTEGER i,j,kg,bi,bj
c     INTEGER iLa
#ifdef LAYERS_THERMODYNAMICS
      INTEGER iTracer
#endif

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C     Loop through 3D layers variables:
        DO kg=1,Nlayers
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx

#ifdef LAYERS_UFLUX
           layers_UH (i,j,kg,bi,bj) = 0. _d 0
#ifdef LAYERS_THICKNESS
           layers_Hw (i,j,kg,bi,bj) = 0. _d 0
           layers_PIw(i,j,kg,bi,bj) = 0. _d 0
           layers_U  (i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_UFLUX */

#ifdef LAYERS_VFLUX
           layers_VH (i,j,kg,bi,bj) = 0. _d 0
#ifdef LAYERS_THICKNESS
           layers_Hs (i,j,kg,bi,bj) = 0. _d 0
           layers_PIs(i,j,kg,bi,bj) = 0. _d 0
           layers_V  (i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THICKNESS */
#endif /* LAYERS_VFLUX */
#ifdef LAYERS_THERMODYNAMICS
           layers_Hc (i,j,kg,bi,bj) = 0. _d 0
           layers_PIc(i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THERMODYNAMICS */
          ENDDO
         ENDDO
        ENDDO
#ifdef LAYERS_THERMODYNAMICS
        DO kg=1,Nlayers-1
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           layers_Hc        (i,j,kg,bi,bj) = 0. _d 0
           layers_TtendSurf (i,j,kg,bi,bj) = 0. _d 0
           layers_TtendDiffh(i,j,kg,bi,bj) = 0. _d 0
           layers_TtendDiffr(i,j,kg,bi,bj) = 0. _d 0
           layers_TtendAdvh (i,j,kg,bi,bj) = 0. _d 0
           layers_TtendAdvr (i,j,kg,bi,bj) = 0. _d 0
           layers_StendSurf (i,j,kg,bi,bj) = 0. _d 0
           layers_StendDiffh(i,j,kg,bi,bj) = 0. _d 0
           layers_StendDiffr(i,j,kg,bi,bj) = 0. _d 0
           layers_StendAdvh(i,j,kg,bi,bj) = 0. _d 0
           layers_StendAdvr(i,j,kg,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
#endif /* LAYERS_THERMODYNAMICS */
#ifdef LAYERS_PRHO_REF
        DO kg=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           prho(i,j,kg,bi,bj) = 0. _d 0
#ifdef LAYERS_THERMODYNAMICS
           layers_alpha(i,j,kg,bi,bj) = 0. _d 0
           layers_beta(i,j,kg,bi,bj) = 0. _d 0
#endif /* LAYERS_THERMODYNAMICS */           
          ENDDO
         ENDDO
        ENDDO
#endif /* LAYERS_PRHO_REF */

#ifdef LAYERS_THERMODYNAMICS
        DO iTracer = 1,2
           DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           layers_surfflux(i,j,1,iTracer,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
         DO kg=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            layers_tottend (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_dfx     (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_dfy     (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_dfr     (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_afx     (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_afy     (i,j,kg,iTracer,bi,bj) = 0. _d 0
            layers_afr     (i,j,kg,iTracer,bi,bj) = 0. _d 0
#ifdef SHORTWAVE_HEATING
            layers_sw      (i,j,kg,1      ,bi,bj) = 0. _d 0
#endif /* SHORTWAVE_HEATING */
           ENDDO
          ENDDO
         ENDDO
        ENDDO
#endif /* LAYERS_THERMODYNAMICS */

C- end bi,bj loops
       ENDDO
      ENDDO

#endif /* ALLOW_LAYERS */

      RETURN
      END
