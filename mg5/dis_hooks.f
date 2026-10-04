C     evgen-benchmark: the DIS generation cuts and the DIS scale for MG5_aMC.
C
C     These two routines REPLACE the upstream stubs of the same names in a
C     generated process directory's SubProcesses/dummy_fct.f; everything else
C     in that file is left as MG5 wrote it.  install_hooks.py does the splice
C     and substitutes @Q2MIN@ / @YMIN@ / @YMAX@ from analysis/selection.py, so
C     the generated sample and the analysis cannot drift apart.
C
C     WHY THESE TWO HOOKS AND NOT THE RUN CARD.
C
C     MG5 has no Q2 or y cut: its run card offers collider cuts (pt, eta,
C     Delta R) only, so the fiducial region has to be imposed in dummy_cuts.
C     And the run card's DEFAULT `dynamical_scale_choice = -1` sets the scale
C     by CLUSTERING the external states; measured on this process it ran at
C     0.787 x Q on average and down to 0.32 x Q, which at these low Q2 cost a
C     FACTOR TWO on the neutral-current cross section with nothing in the
C     output to say so.  `dynamical_scale_choice = 0` routes both muR and muF
C     through user_dynamical_scale below, which returns Q.
C
C     MOMENTUM LABELS.  The process is generated as `mu- q > mu- q` (or
C     `vm qi > mu- qo`), so leg 1 is the incoming lepton, leg 2 the incoming
C     parton, leg 3 the outgoing lepton and leg 4 the outgoing parton, in
C     every subprocess.  Metric (+,-,-,-).
C
C     y NEEDS NO BEAM RECORD: with the incoming parton collinear, p = x P, the
C     x cancels in (p.q)/(p.k), so the parton momentum can be used directly.

      logical FUNCTION dummy_cuts(P)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      REAL*8 P(0:3,nexternal)
      double precision q(0:3), qsq, yy, pq, pk
      integer i
      double precision q2min, ymin, ymax
      parameter (q2min = @Q2MIN@d0)
      parameter (ymin  = @YMIN@d0)
      parameter (ymax  = @YMAX@d0)

      dummy_cuts = .false.

      do i = 0, 3
         q(i) = P(i,1) - P(i,3)
      enddo
      qsq = -(q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3))
      if (qsq .le. q2min) return

      pq = P(0,2)*q(0)   - P(1,2)*q(1)   - P(2,2)*q(2)   - P(3,2)*q(3)
      pk = P(0,2)*P(0,1) - P(1,2)*P(1,1) - P(2,2)*P(2,1) - P(3,2)*P(3,1)
      if (pk .le. 0d0) return
      yy = pq / pk
      if (yy .le. ymin .or. yy .ge. ymax) return

      dummy_cuts = .true.
      return
      end

      double precision function user_dynamical_scale(P)
C     muR = muF = Q, floored at the generation cut so the parton
C     distributions are never asked below the region the benchmark evaluates
C     them in (the standing scale-floor rule).  Reached only with
C     `dynamical_scale_choice = 0` in the run card; choice 0 feeds BOTH
C     scales, so this one routine sets them.
      implicit none
      include 'nexternal.inc'
      double precision P(0:3, nexternal)
      include 'genps.inc'
      include 'vector.inc'
      include 'run.inc'
      double precision q(0:3), qsq
      integer i
      double precision q2min
      parameter (q2min = @Q2MIN@d0)

      do i = 0, 3
         q(i) = P(i,1) - P(i,3)
      enddo
      qsq = -(q(0)*q(0) - q(1)*q(1) - q(2)*q(2) - q(3)*q(3))
      user_dynamical_scale = dsqrt(max(q2min, qsq))
      return
      end
