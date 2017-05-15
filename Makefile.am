################################################################
#                                                              #
#      Automake input file                                     #
#                                                              #
################################################################

# Automake defaults to F77LINK.  Using FC to link will work whether FC is an MPI or non-MPI compiler.
F77LD = $(FC)

PPFCCOMPILE = $(FC) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(FPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(FCFLAGS_f)

.F.o:
	$(F77COMPILE) -c -o $@ $<

.f.o:
	$(PPFCCOMPILE) -c -o $@ $<

FC_MODEXT=mod

# directory with all the *.mod files
MODDIR=$(top_builddir)/model
MODOUTDIR=$(MODDIR)
CLEANFILES=$(MODOUTDIR)/*.$(FC_MODEXT) $(libmfix_src:%.f=%.d) $(postmfix_SOURCES:%.f=%.d) libmfix.a.md5 model/mfix.d
DISTCLEANFILES = build-aux/Makefile.usr
DISTCLEANFILES += benchmarks
DISTCLEANFILES += legacy_tests
DISTCLEANFILES += tests
DISTCLEANFILES += tutorials
DISTCLEANFILES += build-aux/fdep

INCS = $(FC_MODINC)$(top_srcdir)/model/include \
       $(FC_MODINC)$(MODDIR)
AM_FCFLAGS = $(INCS) \
             $(FC_MODOUT)$(MODOUTDIR)
AM_LDFLAGS =

SUFFIXES = .d

.DEFAULT_GOAL := mfix$(EXEEXT)

$(builddir)/read_database.$(OBJEXT): $(srcdir)/model/read_database.f
	$(PPFCCOMPILE) $(FC_DEFINE)BURCAT_THR="'$(srcdir)/model/thermochemical/BURCAT.THR'" -c -o $@ $<

################################################################
#                                                              #
#           mfix executable                                    #
#                                                              #
################################################################
bin_PROGRAMS = mfix$(EXEEXT)

mfix_LDADD = libmfix.a
noinst_LIBRARIES = libmfix.a

# mfix.f needs to be in mfix_SOURCES (not just libmfix_a_SOURCES) for gfortran 4.4
mfix_SOURCES = model/mfix.f
mfix_SOURCES += model/read_database.f
mfix_SOURCES += model/read_namelist.f
libmfix_a_SOURCES = model/ODEPACK.F

################################################################
#                                                              #
#           postmfix executable                                #
#                                                              #
################################################################

bin_PROGRAMS += postmfix$(EXEEXT)

postmfix_LDADD = libmfix.a

# Override $(MODOUTDIR) so post_mfix/functions.mod does not conflict with
# the model/functions.mod.  Both functions.mod files should be the same, but
# compiling both functions_mod.f files in  parallel may cause an error
# if they try to write functions.mod to the same path.
post_mfix/functions_mod.o : MODOUTDIR = $(top_builddir)/post_mfix

CLEANFILES += post_mfix/functions.mod

################################################################
#                                                              #
#      EXTRA_DIST:  Non-source and non-binary files that       #
#                   shoud be included in "make dist" tarball   #
#                                                              #
################################################################

EXTRA_DIST =  Bug_Fixes CHANGES DISCLAIMER NOTICE benchmarks build-aux config configure_mfix doc
EXTRA_DIST += legacy_tests mfix_user_guide.pdf readme regression_tests tests tools tutorials visualization_tools

EXTRA_DIST += model/make_mfix
EXTRA_DIST += model/thermochemical/BURCAT.THR
EXTRA_DIST += model/thermochemical/intro.pdf
EXTRA_DIST += model/thermochemical/readme
distclean-local:
	-if [ "$(srcdir)" != "$(builddir)" ] ; then chmod -R u+w $(DISTCLEANFILES); rm -rf build-aux/test.sh $(DISTCLEANFILES); fi

################################################################
#                                                              #
# Use script "fdep" to find dependencies among source files.   #
# Ignoring the following modules:                              #
#  - omp_lib (OpenMP)                                          #
#  - netcdf                                                    #
#                                                              #
################################################################
DEP_IGN := -i omp_lib -i netcdf -i mpi

.f.d:
	@$(MKDIR_P) $(@D)
	$(top_srcdir)/build-aux/fdep $(DEP_IGN) $(INCS) -m $(MODDIR) -o $(@D) $^

# Use generated module dependency rules
CLEANS := clean distclean mostlyclean
@make_ifnotclean@
$(eval -include $(mfix_SOURCES:%.f=%.d) )
$(eval -include $(libmfix_src:%.f=%.d) )
$(eval -include $(postmfix_SOURCES:%.f=%.d) )
@make_endif@

################################################################
#      TAU support
################################################################
TAU_F90=@TAU_F90@
TAU_MAKEFILE=@TAU_MAKEFILE@
export TAU_F90
export TAU_MAKEFILE

################################################################
#      NETCDF support
################################################################

if NETCDF
AM_FCFLAGS += $(FC_MODINC)$(NETCDF_INCDIR) $(FC_DEFINE)NETCDF
AM_FCFLAGS += $(FC_MODINC)$(NETCDFF_INCDIR) $(FC_DEFINE)NETCDFF
AM_LDFLAGS += -L@NETCDF_LIBDIR@ -l@NETCDF_LIBNAME@ -L@NETCDFF_LIBDIR@ -l@NETCDFF_LIBNAME@
endif

################################################################
#                                                              #
#      MKL support (only with Intel Fortran compiler)          #
#                                                              #
################################################################

if MKL
AM_FCFLAGS += -mkl
AM_LDFLAGS += -mkl
else
libmfix_a_SOURCES += model/BLAS.f
libmfix_a_SOURCES += model/DGTSV.f
endif

################################################################
#                                                              #
#      OpenMP and MPI support                                  #
#                                                              #
#      --enable-smp and --enable-dmp options in configure.ac   #
#                                                              #
################################################################

if OPENMP
AM_FCFLAGS += $(OPENMP_FCFLAGS)
AM_LDFLAGS += $(OPENMP_FCFLAGS)
endif

if MPI
AM_FCFLAGS += $(FC_DEFINE)MPI
AM_LDFLAGS += $(FC_DEFINE)MPI
endif

################################################################
#      libmfix_a_SOURCES                                       #
#                                                              #
#      The majority of source files for mfix and postmfix      #
#      are listed here.                                        #
################################################################

libmfix_a_SOURCES += \
	model/include/boundary_conditions.inc \
	model/include/cartesian_grid_namelist.inc \
	model/include/chem_equations.inc \
	model/include/desnamelist.inc \
	model/include/dmp_batch_control.inc \
	model/include/fun_avg.inc \
	model/include/functions.inc \
	model/include/gas_phase.inc \
	model/include/geometry.inc \
	model/include/initial_conditions.inc \
	model/include/internal_surfaces.inc \
	model/include/legacy.inc \
	model/include/netcdf.inc \
	model/include/numerical_params.inc \
	model/include/output_control.inc \
	model/include/physical_params.inc \
	model/include/point_sources.inc \
	model/include/qmomknamelist.inc \
	model/include/run_control.inc \
	model/include/solids_phase.inc \
	model/include/species.inc \
	model/include/species_indices.inc \
	model/include/tfm_solids.inc \
	model/include/usr_hooks.inc \
	model/include/usrnlst.inc \
	model/include/version.inc \
	model/include/xforms.inc

libmfix_src = \
	model/GhdTheory/adjust_eps_ghd.f \
	model/GhdTheory/bulk_viscosity.f \
	model/GhdTheory/calc_d_ghd.f \
	model/GhdTheory/calc_external_forces.f \
	model/GhdTheory/calc_nflux.f \
	model/GhdTheory/chi_ij_GHD.f \
	model/GhdTheory/cooling_rate.f \
	model/GhdTheory/cooling_rate_tc.f \
	model/GhdTheory/dufour_coeff.f \
	model/GhdTheory/ghd.f \
	model/GhdTheory/ghdmassflux.f \
	model/GhdTheory/ghdtheory_mod.f \
	model/GhdTheory/mass_mobility.f \
	model/GhdTheory/ordinary_diff.f \
	model/GhdTheory/partial_elim_ghd.f \
	model/GhdTheory/pressure.f \
	model/GhdTheory/shear_viscosity.f \
	model/GhdTheory/source_ghd_granular_energy.f \
	model/GhdTheory/thermal_conductivity.f \
	model/GhdTheory/thermal_diffusivity.f \
	model/GhdTheory/thermal_mobility.f \
	model/GhdTheory/transport_coeff_ghd.f \
	model/accum_resid.f \
	model/adjust_a_u_g.f \
	model/adjust_a_u_s.f \
	model/adjust_a_v_g.f \
	model/adjust_a_v_s.f \
	model/adjust_a_w_g.f \
	model/adjust_a_w_s.f \
	model/adjust_eps.f \
	model/adjust_leq.f \
	model/adjust_theta.f \
	model/allocate_arrays.f \
	model/ambm_mod.f \
	model/b_m_p_star.f \
	model/bc_mod.f \
	model/bc_phi.f \
	model/bc_theta.f \
	model/bodyforce_mod.f \
	model/cal_d.f \
	model/calc_cell.f \
	model/calc_coeff.f \
	model/calc_d_e.f \
	model/calc_d_n.f \
	model/calc_d_t.f \
	model/calc_dif_g.f \
	model/calc_dif_s.f \
	model/calc_e.f \
	model/calc_gama.f \
	model/calc_grbdry.f \
	model/calc_h.f \
	model/calc_k_cp.f \
	model/calc_k_g.f \
	model/calc_k_s.f \
	model/calc_mflux.f \
	model/calc_mu_g.f \
	model/calc_mu_s.f \
	model/calc_outflow.f \
	model/calc_p_star.f \
	model/calc_resid.f \
	model/calc_s_ddot_s.f \
	model/calc_trd_g.f \
	model/calc_trd_s.f \
	model/calc_u_friction.f \
	model/calc_vol_fr.f \
	model/cartesian_grid/CG_set_bc0.f \
	model/cartesian_grid/CG_set_outflow.f \
	model/cartesian_grid/CG_source_u_g.f \
	model/cartesian_grid/CG_source_u_s.f \
	model/cartesian_grid/CG_source_v_g.f \
	model/cartesian_grid/CG_source_v_s.f \
	model/cartesian_grid/CG_source_w_g.f \
	model/cartesian_grid/CG_source_w_s.f \
	model/cartesian_grid/allocate_cut_cell_arrays.f \
	model/cartesian_grid/allocate_dummy_cut_cell_arrays.f \
	model/cartesian_grid/calc_vort_out.f \
	model/cartesian_grid/cartesian_grid_init_namelist.f \
	model/cartesian_grid/check_data_cartesian.f \
	model/cartesian_grid/cut_cell_preprocessing.f \
	model/cartesian_grid/cutcell_mod.f \
	model/cartesian_grid/dashboard_mod.f \
	model/cartesian_grid/deallocate_cut_cell_arrays.f \
	model/cartesian_grid/define_quadrics.f \
	model/cartesian_grid/dmp_cartesian.f \
	model/cartesian_grid/eval_usr_fct.f \
	model/cartesian_grid/get_alpha.f \
	model/cartesian_grid/get_connectivity.f \
	model/cartesian_grid/get_cut_cell_flags.f \
	model/cartesian_grid/get_cut_cell_volume_area.f \
	model/cartesian_grid/get_delh.f \
	model/cartesian_grid/get_master.f \
	model/cartesian_grid/get_poly_data.f \
	model/cartesian_grid/get_stl_data.f \
	model/cartesian_grid/polygon_mod.f \
	model/cartesian_grid/progress_bar_mod.f \
	model/cartesian_grid/quadric_mod.f \
	model/cartesian_grid/set_Odxyz.f \
	model/cartesian_grid/stl_mod.f \
	model/cartesian_grid/update_dashboard.f \
	model/cartesian_grid/vtk_mod.f \
	model/cartesian_grid/vtk_out.f \
	model/cartesian_grid/write_progress_bar.f \
	model/cdist_mod.f \
	model/check_ab_m.f \
	model/check_bqend.f \
	model/check_convergence.f \
	model/check_data/check_axis.f \
	model/check_data/check_bc_dem.f \
	model/check_data/check_bc_geometry.f \
	model/check_data/check_bc_inflow.f \
	model/check_data/check_bc_outflow.f \
	model/check_data/check_bc_pic.f \
	model/check_data/check_bc_walls.f \
	model/check_data/check_boundary_conditions.f \
	model/check_data/check_chemical_rxns.f \
	model/check_data/check_dmp_prereqs.f \
	model/check_data/check_gas_phase.f \
	model/check_data/check_geometry.f \
	model/check_data/check_geometry_prereqs.f \
	model/check_data/check_ic_common_discrete.f \
	model/check_data/check_ic_dem.f \
	model/check_data/check_ic_mppic.f \
	model/check_data/check_initial_conditions.f \
	model/check_data/check_internal_surfaces.f \
	model/check_data/check_numerics.f \
	model/check_data/check_odepack_stiff_chem.f \
	model/check_data/check_output_control.f \
	model/check_data/check_point_sources.f \
	model/check_data/check_run_control.f \
	model/check_data/check_solids_common_all.f \
	model/check_data/check_solids_common_discrete.f \
	model/check_data/check_solids_continuum.f \
	model/check_data/check_solids_dem.f \
	model/check_data/check_solids_model_prereqs.f \
	model/check_data/check_solids_mppic.f \
	model/check_data/check_solids_phases.f \
	model/check_data_20.f \
	model/check_data_30.f \
	model/check_mod.f \
	model/check_plane.f \
	model/chem/stiff_chem_dbg_mod.f \
	model/chem/stiff_chem_maps_mod.f \
	model/chem/stiff_chem_mod.f \
	model/chem/stiff_chem_rrates.f \
	model/chem/stiff_chem_stats_mod.f \
	model/chischeme_mod.f \
	model/cn_extrapol.f \
	model/coeff_mod.f \
	model/constant_mod.f \
	model/conv_dif_phi.f \
	model/conv_dif_u_g.f \
	model/conv_dif_u_s.f \
	model/conv_dif_v_g.f \
	model/conv_dif_v_s.f \
	model/conv_dif_w_g.f \
	model/conv_dif_w_s.f \
	model/conv_pp_g.f \
	model/conv_rop.f \
	model/conv_rop_g.f \
	model/conv_rop_s.f \
	model/conv_source_epp.f \
	model/copy_a.f \
	model/corner.f \
	model/corner_mod.f \
	model/correct_0.f \
	model/correct_1.f \
	model/dbg_mod.f \
	model/deprecated.f \
	model/des/calc_collision_wall_mod.f \
  model/des/calc_des_2fluid.f \
	model/des/calc_drag_des.f \
	model/des/calc_epg_des.f \
	model/des/calc_force_dem.f \
  model/des/calc_gamma_des.f \
	model/des/calc_grad_des.f \
	model/des/calc_interp_weights.f \
	model/des/calc_pg_grad.f \
	model/des/calc_rrates_des.f \
	model/des/calc_thermo_des.f \
	model/des/cfassign.f \
	model/des/cffctowall.f \
	model/des/cfnewvalues.f \
	model/des/cfrelvel.f \
	model/des/cfslide.f \
	model/des/cfupdateold.f \
	model/des/cfwallcontact.f \
	model/des/cfwallposvel.f \
	model/des/check_cell_movement.f \
	model/des/comp_mean_fields.f \
	model/des/comp_mean_fields0.f \
	model/des/comp_mean_fields1.f \
	model/des/derived_types_mod.f \
	model/des/des_allocate_mod.f \
	model/des/des_bc_mod.f \
	model/des/des_cluster_mod.f \
	model/des/des_functions.f \
	model/des/des_granular_temperature.f \
	model/des/des_init_arrays.f \
	model/des/des_init_namelist.f \
	model/des/des_physical_prop.f \
	model/des/des_reaction_model.f \
	model/des/des_rxns_mod.f \
	model/des/des_thermo_cond_mod.f \
	model/des/des_thermo_conv.f \
	model/des/des_thermo_mod.f \
	model/des/des_thermo_newvalues.f \
	model/des/des_thermo_rad.f \
	model/des/des_time_march.f \
	model/des/desgrid_mod.f \
	model/des/desmpi_mod.f \
	model/des/desmpi_wrapper_mod.f \
	model/des/dif_phi_bc_des.f \
	model/des/dif_phi_des.f \
	model/des/dif_phi_source_des.f \
	model/des/diffuse_mean_fields.f \
	model/des/discretelement_mod.f \
	model/des/drag_gp_des.f \
	model/des/drag_gs_des0.f \
	model/des/drag_gs_des1.f \
	model/des/drag_ss_dem_noninterp.f \
	model/des/gas_drag.f \
	model/des/generate_particles_mod.f \
	model/des/init_settling_dem.f \
	model/des/interpolation_mod.f \
	model/des/layout_mi_dem.f \
	model/des/list.f \
	model/des/make_arrays_des.f \
	model/des/mass_inflow_dem.f \
	model/des/mass_outflow_dem.f \
	model/des/mpi_comm_des_mod.f \
	model/des/mpi_funs_des_mod.f \
	model/des/mpi_init_des_mod.f \
	model/des/mpi_node_des_mod.f \
	model/des/mpi_pack_des_mod.f \
	model/des/mpi_unpack_des_mod.f \
	model/des/multisap.f \
	model/des/neighbour.f \
	model/des/nsquare.f \
	model/des/pair_manager.f \
	model/des/particle_filter_mod.f \
	model/des/particles_in_cell.f \
	model/des/pic/apply_wall_bc_pic.f \
	model/des/pic/calc_dtpic.f \
	model/des/pic/calc_ps_grad_pic.f \
	model/des/pic/calc_ps_pic.f \
	model/des/pic/integrate_time_pic.f \
	model/des/pic/interpolate_pic.f \
	model/des/pic/mass_inflow_pic.f \
	model/des/pic/mass_outflow_pic.f \
	model/des/pic/mfix_pic_mod.f \
	model/des/pic/pic_bc_mod.f \
	model/des/pic/pic_routines.f \
	model/des/pic/report_stats_pic.f \
	model/des/pic/set_bc_pic.f \
	model/des/pic/set_bc_pic_mi.f \
	model/des/pic/set_bc_pic_mo.f \
	model/des/pic/time_march_pic.f \
	model/des/randomno_mod.f \
	model/des/read_particle_input.f \
	model/des/read_res0_des.f \
	model/des/read_res1_des_mod.f \
	model/des/resize.f \
  model/des/rxns_gs_des1.f \
	model/des/sap.f \
	model/des/sendrecvnode_mod.f \
	model/des/set_bc_dem.f \
	model/des/set_bc_dem_mi.f \
	model/des/set_bc_dem_mo.f \
	model/des/set_filter_des.f \
	model/des/set_geometry_des.f \
	model/des/set_ic_dem.f \
	model/des/set_phase_index.f \
	model/des/solid_drag.f \
	model/des/stl_dbg_des_mod.f \
	model/des/stl_functions_des_mod.f \
	model/des/stl_preproc_des_mod.f \
	model/des/usr0_des.f \
	model/des/usr1_des.f \
	model/des/usr2_des.f \
	model/des/usr3_des.f \
	model/des/usr4_des.f \
	model/des/usr_rates_des.f \
	model/des/vtp_mod.f \
	model/des/write_des_data.f \
	model/des/write_res0_des.f \
	model/des/write_res1_des_mod.f \
	model/dgtsl.f \
	model/dif_u_is.f \
	model/dif_v_is.f \
	model/dif_w_is.f \
	model/discretization_mod.f \
	model/display_resid.f \
	model/dmp_modules/compar_mod.f \
	model/dmp_modules/dbg_util_mod.f \
	model/dmp_modules/debug_mod.f \
	model/dmp_modules/gridmap_mod.f \
	model/dmp_modules/mpi_utility_mod.f \
	model/dmp_modules/parallel_mpi_mod.f \
	model/dmp_modules/sendrecv3_mod.f \
	model/dmp_modules/sendrecv_mod.f \
	model/dqmom/gaussj.f \
	model/dqmom/odeint.f \
	model/dqmom/rkck.f \
	model/dqmom/rkqs.f \
	model/dqmom/source_population_eq.f \
	model/dqmom/usr_dqmom.f \
	model/drag_gs.f \
	model/drag_mod.f \
	model/drag_ss.f \
	model/energy_mod.f \
	model/eos_mod.f \
	model/error_manager_mod.f \
	model/exchange.f \
	model/exit.f \
	model/fldvar_mod.f \
	model/flow_to_vel.f \
	model/fun_avg_mod.f \
	model/function3_mod.f \
	model/functions_mod.f \
	model/funits_mod.f \
	model/geometry_mod.f \
	model/get_bc_area.f \
	model/get_data.f \
	model/get_eq.f \
	model/get_flow_bc.f \
	model/get_hloss.f \
	model/get_is.f \
	model/get_philoss.f \
	model/get_ps.f \
	model/get_smass.f \
	model/get_stats.f \
	model/get_walls_bc.f \
	model/ic_mod.f \
	model/in_binary_512_mod.f \
	model/in_binary_512i_mod.f \
	model/indices_mod.f \
	model/init_ab_m.f \
	model/init_fvars.f \
	model/init_namelist.f \
	model/init_resid.f \
	model/is_mod.f \
	model/iterate.f \
	model/k_epsilon_prop.f \
	model/kintheory_mod.f \
	model/kintheory_u_s.f \
	model/kintheory_v_s.f \
	model/kintheory_w_s.f \
	model/leq_bicgs.f \
	model/leq_bicgst.f \
	model/leq_cg.f \
	model/leq_gmres.f \
	model/leq_sor.f \
	model/leqsol_mod.f \
	model/location_check.f \
	model/machine_mod.f \
	model/main.f \
	model/make_upper_case.f \
	model/mark_phase_4_cor.f \
	model/mfix_netcdf_mod.f \
	model/mflux_mod.f \
	model/mms_mod.f \
	model/mod_bc_i.f \
	model/mod_bc_j.f \
	model/mod_bc_k.f \
	model/open_file.f \
	model/open_files.f \
	model/out_array.f \
	model/out_array_c.f \
	model/out_array_k.f \
	model/out_array_kc.f \
	model/out_bin_512.f \
	model/out_bin_512i.f \
	model/out_bin_512r.f \
	model/out_bin_r.f \
	model/output_manager.f \
	model/output_mod.f \
	model/parallel_mod.f \
	model/param1_mod.f \
	model/param_mod.f \
	model/parse_line.f \
	model/parse_mod.f \
	model/parse_resid_string.f \
	model/parse_rxn.f \
	model/partial_elim.f \
	model/pgcor_mod.f \
	model/physical_prop.f \
	model/physprop_mod.f \
	model/ps_mod.f \
	model/pscor_mod.f \
	model/qmomk/qmom_kinetic_equation_mod.f \
	model/qmomk/qmomk_allocate_arrays.f \
	model/qmomk/qmomk_bc_mod.f \
	model/qmomk/qmomk_collision_mod.f \
	model/qmomk/qmomk_fluxes_mod.f \
	model/qmomk/qmomk_gas_drag.f \
	model/qmomk/qmomk_init_bc.f \
	model/qmomk/qmomk_init_namelist.f \
	model/qmomk/qmomk_initial_conditions.f \
	model/qmomk/qmomk_make_arrays.f \
	model/qmomk/qmomk_parameters_mod.f \
	model/qmomk/qmomk_quadrature_mod.f \
	model/qmomk/qmomk_read_restart.f \
	model/qmomk/qmomk_set_bc.f \
	model/qmomk/qmomk_time_march.f \
	model/qmomk/qmomk_tools_mod.f \
	model/qmomk/qmomk_write_restart.f \
	model/radial_vel_correction.f \
	model/rdf_mod.f \
	model/read_res1.f \
	model/reinitialize.f \
	model/remove_comment.f \
	model/reset_new.f \
	model/residual_mod.f \
	model/rrates.f \
	model/rrates0.f \
	model/rrates_init.f \
	model/run_mod.f \
	model/rxn_com_mod.f \
	model/rxns_mod.f \
	model/scalar_prop.f \
	model/scalars_mod.f \
	model/scales_mod.f \
	model/set_bc0.f \
	model/set_bc1.f \
	model/set_bc_flow.f \
	model/set_constprop.f \
	model/set_flags.f \
	model/set_fluidbed_p.f \
	model/set_geometry.f \
	model/set_geometry1.f \
	model/set_ic.f \
	model/set_icbc_flags.f \
	model/set_increments.f \
	model/set_increments3.f \
	model/set_index1.f \
	model/set_index1a.f \
	model/set_index1a3.f \
	model/set_l_scale.f \
	model/set_max2.f \
	model/set_mw_mix_g.f \
	model/set_outflow.f \
	model/set_param.f \
	model/set_ps.f \
	model/set_ro_g.f \
	model/set_ro_s.f \
	model/set_wall_bc.f \
	model/shift_dxyz.f \
	model/solids_pressure_mod.f \
	model/solve_continuity.f \
	model/solve_energy_eq.f \
	model/solve_epp.f \
	model/solve_granular_energy.f \
	model/solve_k_epsilon_eq.f \
	model/solve_lin_eq.f \
	model/solve_pp_g.f \
	model/solve_scalar_eq.f \
	model/solve_species_eq.f \
	model/solve_vel_star.f \
	model/source_granular_energy.f \
	model/source_phi.f \
	model/source_pp_g.f \
	model/source_rop_g.f \
	model/source_rop_s.f \
	model/source_u_g.f \
	model/source_u_s.f \
	model/source_v_g.f \
	model/source_v_s.f \
	model/source_w_g.f \
	model/source_w_s.f \
	model/time_step.f \
	model/tau_g_mod.f \
	model/tau_s_mod.f \
	model/tau_u_g.f \
	model/tau_u_s.f \
	model/tau_v_g.f \
	model/tau_v_s.f \
	model/tau_w_g.f \
	model/tau_w_s.f \
	model/test_lin_eq.f \
	model/thermochemical/get_values.f \
	model/thermochemical/read_thermochemical_mod.f \
	model/time_cpu_mod.f \
	model/toleranc_mod.f \
	model/trace_mod.f \
	model/transport_prop.f \
	model/turb_mod.f \
	model/undef_2_0.f \
	model/under_relax.f \
	model/update_old.f \
	model/ur_facs_mod.f \
	model/usr0.f \
	model/usr1.f \
	model/usr2.f \
	model/usr3.f \
	model/usr_drag.f \
	model/usr_init_namelist.f \
	model/usr_mod.f \
	model/usr_prop_mod.f \
	model/usr_properties.f \
	model/usr_rates.f \
	model/usr_sources.f\
	model/usr_src_mod.f \
	model/usr_write_out0.f \
	model/usr_write_out1.f \
	model/utilities_mod.f \
	model/vavg_u_g.f \
	model/vavg_u_s.f \
	model/vavg_v_g.f \
	model/vavg_v_s.f \
	model/vavg_w_g.f \
	model/vavg_w_s.f \
	model/vf_gs_x.f \
	model/vf_gs_y.f \
	model/vf_gs_z.f \
	model/visc_g_mod.f \
	model/visc_s_mod.f \
	model/vshear_mod.f \
	model/vtc_scalar.f \
	model/write_ab_m.f \
	model/write_ab_m_var.f \
	model/write_error.f \
	model/write_header.f \
	model/write_out0.f \
	model/write_out1.f \
	model/write_out3.f \
	model/write_res0.f \
	model/write_res1.f \
	model/write_spx0.f \
	model/write_spx1.f \
	model/write_table.f \
	model/write_usr0.f \
	model/write_usr1.f \
	model/xerbla.f \
	model/xsi_mod.f \
	model/zero_norm_vel.f

libmfix_a_SOURCES += $(libmfix_src)

################################################################
#                                                              #
#      postmfix                                                #
#                                                              #
################################################################

postmfix_SOURCES = \
	model/read_database.f \
	model/read_namelist.f \
	post_mfix/any_more_data.f \
	post_mfix/calc_cell2.f \
	post_mfix/calc_corr_01.f \
	post_mfix/calc_corr_type_1.f \
	post_mfix/calc_distance.f \
	post_mfix/calc_dpody.f \
	post_mfix/calc_ep_g.f \
	post_mfix/calc_quantities.f \
	post_mfix/calc_ro_g.f \
	post_mfix/calc_vol.f \
	post_mfix/correl_mod.f \
	post_mfix/deallocate_arrays.f \
	post_mfix/examine_data.f \
	post_mfix/f_init_data.f \
	post_mfix/file_handle.f \
	post_mfix/finit.f \
	post_mfix/flow_gx.f \
	post_mfix/flow_gy.f \
	post_mfix/flow_gz.f \
	post_mfix/flow_sx.f \
	post_mfix/flow_sy.f \
	post_mfix/flow_sz.f \
	post_mfix/functions_mod.f \
	post_mfix/gas_flux.f \
	post_mfix/get_file_name.f \
	post_mfix/get_file_status.f \
	post_mfix/get_index.f \
	post_mfix/get_location.f \
	post_mfix/get_mu_s.f \
	post_mfix/get_same_time.f \
	post_mfix/get_selection.f \
	post_mfix/get_substr.f \
	post_mfix/granular_qty.f \
	post_mfix/header_main.f \
	post_mfix/ik_avg.f \
	post_mfix/ik_avg_out.f \
	post_mfix/interp_res.f \
	post_mfix/main_f.f \
	post_mfix/ornl_corr.f \
	post_mfix/ornl_corr_c.f \
	post_mfix/ornl_filt.f \
	post_mfix/ornl_filt_c.f \
	post_mfix/ornl_ft.f \
	post_mfix/ornl_ft_c.f \
	post_mfix/ornl_header.f \
	post_mfix/ornl_pca.f \
	post_mfix/ornl_stats.f \
	post_mfix/ornl_stats_c.f \
	post_mfix/ornl_sym.f \
	post_mfix/ornl_util.f \
	post_mfix/ornl_zone.f \
	post_mfix/out_from_res.f \
	post_mfix/out_from_spx.f \
	post_mfix/out_spec_time.f \
	post_mfix/out_time.f \
	post_mfix/paralleldata_mod.f \
	post_mfix/post3d_mod.f \
	post_mfix/post_precision_mod.f \
	post_mfix/print_out.f \
	post_mfix/read_res0.f \
	post_mfix/read_spx0.f \
	post_mfix/read_spx1.f \
	post_mfix/res_from_spx.f \
	post_mfix/seek_time.f \
	post_mfix/select_spx_rec.f \
	post_mfix/set_dollar.f \
	post_mfix/set_read_spx.f \
	post_mfix/sol_flux.f \
	post_mfix/strcmp.f \
	post_mfix/streqs.f \
	post_mfix/time_avg.f \
	post_mfix/usr_init_namelist.f \
	post_mfix/usr_input.f \
	post_mfix/usr_input_mod.f \
	post_mfix/usr_post.f \
	post_mfix/usr_write_out1.f

cases := $(shell find . -name mfix.dat)

allcases:
	for case in $(cases); do \
	  echo "BUILDING: $$case"; \
	  cp build-aux/Makefile.usr $$(dirname $$case)/Makefile || exit; \
	  $(MAKE) -C $$(dirname $$case) mfix$(EXEEXT) || exit; \
	done

TESTS = build-aux/test.sh
