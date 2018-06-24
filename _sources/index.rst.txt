.. pytheas documentation master file, created by
   sphinx-quickstart on Tue Jun  6 20:04:57 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.




.. title:: Home


   
Pytheas is a python package for creating, running and postprocessing
electrodynamic simulations. It is based on the open source onelab 
bundle for creating geometries, meshing (gmsh) and solving 
the underlying partial differential equations with the finite 
element method (getdp).
It features built in models of:

- periodic media in 2D and 3D with computation of diffraction efficiencies
- scattering analysis in 2D and 3D
- bloch mode analysis of metamaterials
- treatment of open geometries with perfectly matched layers
- tools to define arbitrary permittivity distributions
- quasi-normal mode analysis
- two scale convergence homogenization
- tools for topology optimization in 2D
- built-in refractive index database


.. raw:: html
  
  <div class="clearer">
    <div class="row">
      <div class="col-md-6">
        <p class="button">
          <a type="button" class="btn btn-primary" href="reference.html">
            Read the doc
          </a>
          </p>
          </div> 
          <div class="col-md-6">
          <p class="button">
          <a type="button" class="btn btn-primary" href="auto_examples/index.html">
            Browse examples
          </a>
          </p>
        </div> 
    </div>
  </div>



.. toctree::
   :maxdepth: 2
   :hidden:

   reference
   auto_examples/index



.. Indices and tables
.. ==================
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
