# script to generate all the figures

exec(open('dispersion_relation_wake.py').read())
exec(open('ship_wake.py').read())
exec(open('wake_according_to_Raphael.py').read())
exec(open('plot_sillage_experimental.py').read())
exec(open('wavelength.py').read())


import matlab.engine

eng = matlab.engine.start_matlab()
eng.field_modal(nargout=0)
eng.kreal_modal(nargout=0)
eng.kcomplex_field(nargout=0)
eng.kcomplex_comparison(nargout=0)
eng.comparacion_modal_bc_efectiva(nargout=0)
eng.comparacion_reld_sp_eigenvalues_adim(nargout=0)
eng.quit()

print('Done.')

