#include "GMVIO_cell.h"

// anonymous namespace to hold local data
namespace
{
  /**
   * Defines mapping from libMesh element types to GMV element types.
   */
  struct elementDefinition {
    std::string label;
    std::vector<unsigned int> nodes;
  };


  // maps from a libMesh element type to the proper
  // GMV elementDefinition.  Placing the data structure
  // here in this anonymous namespace gives us the
  // benefits of a global variable without the nasty
  // side-effects
  std::map<ElemType, elementDefinition> eletypes;



  // ------------------------------------------------------------
  // helper function to initialize the eletypes map
  void init_eletypes ()
  {
    if (eletypes.empty())
      {
	// This should happen only once.  The first time this method
	// is called the eletypes data struture will be empty, and
	// we will fill it.  Any subsequent calls will find an initialized
	// eletypes map and will do nothing.

	//==============================
	// setup the element definitions
	elementDefinition eledef;

	// use "swap trick" from Scott Meyer's "Effective STL" to initialize
	// eledef.nodes vector
	
	// EDGE2
	{
	  eledef.label = "line 2";
	  const unsigned int nodes[] = {0,1};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	
	  eletypes[EDGE2] = eledef;
	}
  
	// LINE3
	{
	  eledef.label = "line 3";
	  const unsigned int nodes[] = {0,1,2};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[EDGE3] = eledef;
	}
      
	// TRI3
	{
	  eledef.label = "tri3 3";
	  const unsigned int nodes[] = {0,1,2};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[TRI3] = eledef;
	}
      
	// TRI6
	{
	  eledef.label = "6tri 6";
	  const unsigned int nodes[] = {0,1,2,3,4,5};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[TRI6] = eledef;
	}
      
	// QUAD4
	{
	  eledef.label = "quad 4";
	  const unsigned int nodes[] = {0,1,2,3};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[QUAD4] = eledef;
	}
      
	// QUAD8, QUAD9
	{
	  eledef.label = "8quad 8";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[QUAD8] = eledef;
	  eletypes[QUAD9] = eledef;
	}
      
	// HEX8
	{
	  eledef.label = "phex8 8";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[HEX8] = eledef;
	}
      
	// HEX20, HEX27
	{
	  eledef.label = "phex20 20";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[HEX20] = eledef;
	  eletypes[HEX27] = eledef;
	}
      
	// TET4
	{
	  eledef.label = "tet 4";
	  const unsigned int nodes[] = {0,1,2,3};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);

	  eletypes[TET4] = eledef;
	}
      
	// TET10
	{
	  eledef.label = "tet10 10";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,9};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[TET10] = eledef;
	}
      
	// PRISM6
	{
	  eledef.label = "pprism6 6";
	  const unsigned int nodes[] = {0,1,2,3,4,5};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[PRISM6] = eledef;
	}
      
	// PRISM15, PRISM18
	{
	  eledef.label = "pprism15 15";
	  const unsigned int nodes[] = {0,1,2,3,4,5,6,7,8,12,13,14,9,10,11};
	  const unsigned int nnodes = sizeof(nodes)/sizeof(nodes[0]);
	  std::vector<unsigned int>(nodes, nodes+nnodes).swap(eledef.nodes);
	  
	  eletypes[PRISM15] = eledef;
	  eletypes[PRISM18] = eledef;
	}
	//==============================      
      }
  }
  
} // end anonymous namespace


void GMVIO_cell::write_ascii_cell_data (const std::string& fname,
				  const std::vector<Number>& v1,
				  const std::vector<std::string>& solution_names1)

{
  const std::vector<Number>* v = &v1;
  const std::vector<std::string>* solution_names = &solution_names1;
  
 
  // Open the output file stream
  std::ofstream out (fname.c_str());
 
  assert (out.good());
  
 
  // Get a reference to the mesh
  const MeshBase& mesh = MeshOutput<MeshBase>::mesh();


  // Begin interfacing with the GMV data file
  {
    out << "gmvinput ascii\n\n";

    // write the nodes
    out << "nodes " << mesh.n_nodes() << "\n";
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(0) << " ";
    out << "\n";
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(1) << " ";
    out << "\n";
    
    for (unsigned int v=0; v<mesh.n_nodes(); v++)
      out << mesh.point(v)(2) << " ";
    out << "\n\n";
  }

  {
    // write the connectivity
    out << "cells " << mesh.n_active_elem() << "\n";
    
    // initialize the eletypes map
    init_eletypes();

    MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end = mesh.active_elements_end(); 
    
    for ( ; it != end; ++it)
      {
        const Elem* elem = *it;

	// Make sure we have a valid entry for
	// the current element type.
	assert (eletypes.count(elem->type()));

        const elementDefinition& ele = eletypes[elem->type()];
	
        out << ele.label << "\n";
        for (unsigned int i=0; i < ele.nodes.size(); i++)
          out << elem->node(ele.nodes[i])+1 << " ";
        out << "\n";
      }
    out << "\n";
  }
  
  // optionally write the partition information
  if (this->partitioning())
    {
      out << "material "
          << mesh.n_partitions()
          << " 0\n";

      for (unsigned int proc=0; proc<mesh.n_partitions(); proc++)
        out << "proc_" << proc << "\n";
      
      MeshBase::const_element_iterator       it  = mesh.active_elements_begin();
      const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

      for ( ; it != end; ++it)
        out << (*it)->processor_id()+1 << "\n";
      out << "\n";
    }





  // NEW :  write ELEMENT data !!!!


  // optionally write the data
  if ((solution_names != NULL) && (v != NULL))
    {      
      const unsigned int n_vars = solution_names->size();

      if (!(v->size() == mesh.n_active_elem()*n_vars))
        std::cerr << "ERROR: v->size()=" << v->size()
                  << ", mesh.n_active_elem()=" << mesh.n_active_elem()
                  << ", n_vars=" << n_vars
                  << ", mesh.n_active_elem()*n_vars=" << mesh.n_active_elem()*n_vars
                  << "\n";
      
      assert (v->size() == mesh.n_active_elem()*n_vars);

      out << "variable" << "\n";

      for (unsigned int c=0; c<n_vars; c++)
        {

#ifdef USE_COMPLEX_NUMBERS

          // in case of complex data, write _three_ data sets
          // for each component

          // this is the real part
          out << "r_" << (*solution_names)[c] << " 0\n";
	  
          for (unsigned int n=0; n<mesh.n_active_elem(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c].real() << " ";

          out << "\n\n";

          // this is the imaginary part
          out << "i_" << (*solution_names)[c] << " 0\n";
	  
          for (unsigned int n=0; n<mesh.n_active_elem(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c].imag() << " ";

          out << "\n\n";

          // this is the magnitude
          out << "a_" << (*solution_names)[c] << " 0\n";
          for (unsigned int n=0; n<mesh.n_active_elem(); n++)
            out << std::setprecision(10)
                << std::abs((*v)[n*n_vars + c]) << " ";

          out << "\n\n";

#else

          out << (*solution_names)[c] << " 0\n"; //  0  for  cell (element) data  !!!!!!
	  

	  const unsigned int n_active_elements = mesh.n_active_elem();

          for (unsigned int n=0; n < n_active_elements ; n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c] << " ";
/*	  
          for (unsigned int n=0; n<mesh.n_nodes(); n++)
            out << std::setprecision(10) << (*v)[n*n_vars + c] << " ";

*/
          out << "\n\n";

#endif
        }
      
      out << "endvars\n";
    }
  out << "\n endgmv \n";
}
