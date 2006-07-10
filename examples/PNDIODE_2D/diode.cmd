Title "Untitled"

Definitions {

  # Refinement regions 
    Refinement "Default Region"
    {
      MaxElementSize = (10 10 )
      MinElementSize = (5 1 )
    }
    Refinement "NoName_4"
    {
      MaxElementSize = (1 1 )
      MinElementSize = (0.5 0.005 )
    }
    Refinement "NoName_8"
    {
      MaxElementSize = (10 0.1 )
      MinElementSize = (0.02 0.001 )
    }
    Refinement "NoName_3"
    {
      MaxElementSize = (1 1 )
      MinElementSize = (0.5 0.005 )
    }
    Refinement "NoName_6"
    {
      MaxElementSize = (25 0.01 )
      MinElementSize = (0.5 0.001 )
    }
    Refinement "NoName_5"
    {
      MaxElementSize = (0.5 0.1 )
      MinElementSize = (0.02 0.005 )
    }
    Refinement "NoName_7"
    {
      MaxElementSize = (0.1 0.05 )
      MinElementSize = (0.02 0.005 )
    }
    Refinement "NoName_2"
    {
      MaxElementSize = (5 2.5 )
      MinElementSize = (1 0.005 )
    }

  # Profiles 
    Constant "n_dop"
    {
      Species = "ArsenicActiveConcentration"
      Value = 1e+18
    }
    Constant "p_dop"
    {
      Species = "BoronActiveConcentration"
      Value = 1e+18
    }
}

Placements {

  # Refinement regions 
    Refinement "Default Region"
    {
      Reference = "Default Region"
      # Default region 
    }
    Refinement "NoName_2"
    {
      Reference = "NoName_2"
      RefineWindow = rectangle [( -2.7162 35 ) , ( 102.72 77.47 )]
    }
    Refinement "NoName_8"
    {
      Reference = "NoName_8"
      RefineWindow = rectangle [( -0.76923 45 ) , ( 100.77 55 )]
    }
    Refinement "NoName_3"
    {
      Reference = "NoName_3"
      RefineWindow = rectangle [( 69 47.97 ) , ( 85 50.884 )]
    }
    Refinement "NoName_5"
    {
      Reference = "NoName_5"
      RefineWindow = rectangle [( 58 49 ) , ( 70 51 )]
    }
    Refinement "NoName_7"
    {
      Reference = "NoName_7"
      RefineWindow = rectangle [( 59.5 49.2 ) , ( 62 50 )]
    }
    Refinement "NoName_6"
    {
      Reference = "NoName_6"
      RefineWindow = rectangle [( -0.7132 49.9 ) , ( 61.381 50.1 )]
    }
    Refinement "NoName_4"
    {
      Reference = "NoName_4"
      RefineWindow = rectangle [( 30 72.917 ) , ( 45 75.831 )]
    }

  # Profiles 
    Constant "n_dop"
    {
      Reference = "n_dop"
      EvaluateWindow
      {
        Element = rectangle [( 0 0 ) , ( 100 50 )]
        DecayLength = 0
      }
    }

    Constant "p_dop"
    {
      Reference = "p_dop"
      EvaluateWindow
      {
        Element = rectangle [( 0 50 ) , ( 60 75 )]
        DecayLength = 0
      }
    }

}

Offsetting {

#steps to perform
usebox = 0

#global settings:
noffset {
hloc= 0
factor= 1.3
subdivide= 0
}
noffset {
maxedgelength= 1e+30
terminateline= 3
maxlevel= 200
}
boundary {
hglob= 1e+08
}

#isolines:
}
