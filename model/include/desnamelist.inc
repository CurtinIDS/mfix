       NAMELIST / DES_INPUT_DATA /&

! USE discretelement
        DES_ONEWAY_COUPLED, DES_REPORT_MASS_INTERP, DES_INTG_METHOD, &
        PARTICLES, NFACTOR, USE_VDH_DEM_MODEL, &
        GENER_PART_CONFIG, DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO, DES_USR_VAR_SIZE,&

        DES_INTERP_SCHEME, DES_INTERP_WIDTH, DES_INTERP_ON, &
        DES_INTERP_MEAN_FIELDS, DES_EXPLICITLY_COUPLED, &
        DES_DIFFUSE_WIDTH, &

! Used for grid/cell based search methods
        DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, DESGRIDSEARCH_KMAX, &

! Particle properties
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, &
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!
! Additional properties for Hertzian model:
        DES_COLL_MODEL, EW_YOUNG, VW_POISSON, E_YOUNG, V_POISSON,&

! Actual material properties for particles/wall.  Used for computing
! conduction correction terms to account for artificial softening
  	EW_YOUNG_ACTUAL, E_YOUNG_ACTUAL, V_POISSON_ACTUAL, &
	VW_POISSON_ACTUAL, &

! Lees Edwards BC (for DEM only)
        DES_LE_BC, DES_LE_REL_VEL, DES_LE_SHEAR_DIR,&


! DES wall BC: wall translational velocity
! These need to be improved to work for arbitrary geoms
        DES_BC_Uw_s, DES_BC_Vw_s, DES_BC_Ww_s,&

! DES flow boundary condition information: J.Musser
        FORCE_ORD_BC, &

! DES output control information
        DEBUG_DES, FOCUS_PARTICLE, VTP_DIR, &
        PRINT_DES_DATA, DES_OUTPUT_TYPE,  &

! DES Cohesion Information: M.Weber
        USE_COHESION,VAN_DER_WAALS, HAMAKER_CONSTANT, &
        WALL_HAMAKER_CONSTANT, &
        VDW_INNER_CUTOFF, VDW_OUTER_CUTOFF, ASPERITIES, &
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,&

! DES Energy Equations: J.Musser
        DES_CONV_CORR, FLPC, DES_Em, DES_MIN_COND_DIST, &

! MP-PIC related flags
         MPPIC_SOLID_STRESS_SNIDER, PSFAC_FRIC_PIC,&
         FRIC_EXP_PIC, FRIC_NON_SING_FAC,&
         CFL_PIC, MPPIC_COEFF_EN1, MPPIC_COEFF_EN2, MPPIC_PDRAG_IMPLICIT, &
         MPPIC_GRAV_TREATMENT, MPPIC_COEFF_EN_WALL, MPPIC_COEFF_ET_WALL, &
         PIC_REPORT_SEEDING_STATS, PIC_REPORT_DELETION_STATS, PIC_REPORT_MIN_EPG, &

! Flag to turn on/off optimizing the list of facets at each des grid cell
         MINIMIZE_DES_FACET_LIST
