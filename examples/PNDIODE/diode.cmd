Title "Untitled"

Definitions {

  # Refinement regions 
    Refinement "Default Region"
    {
      MaxElementSize = (0.2 0.2 )
      MinElementSize = (0.001 0.001 )
      RefineFunction = MaxTransDiff(Variable = "DopingConcentration", Value = 1)
    }
    Refinement "NoName_1"
    {
      MaxElementSize = (0.05 0.2 )
      MinElementSize = (0.001 0.001 )
    }
    Refinement "NoName_2"
    {
      MaxElementSize = (0.05 0.2 )
      MinElementSize = (0.001 0.001 )
    }
    Refinement "NoName_3"
    {
      MaxElementSize = (0.05 0.05 )
      MinElementSize = (0.001 0.001 )
    }
    Refinement "NoName_4"
    {
      MaxElementSize = (0.05 0.05 )
      MinElementSize = (0.001 0.001 )
    }

  # Profiles 
    Constant "ndop"
    {
      Species = "ArsenicActiveConcentration"
      Value = 1e+18
    }
    Constant "pdop"
    {
      Species = "BoronActiveConcentration"
      Value = 1e+19
    }
}

Placements {

  # Refinement regions 
    Refinement "Default Region"
    {
      Reference = "Default Region"
      # Default region 
    }
    Refinement "NoName_1"
    {
      Reference = "NoName_1"
      RefineWindow = rectangle [( -0.89531 0.88516 ) , ( -0.56719 2.0664 )]
    }
    Refinement "NoName_2"
    {
      Reference = "NoName_2"
      RefineWindow = rectangle [( 0.57656 0.87578 ) , ( 0.90938 2.0523 )]
    }
    Refinement "NoName_3"
    {
      Reference = "NoName_3"
      RefineWindow = rectangle [( -0.89531 0.88047 ) , ( 0.9 1.1242 )]
    }
    Refinement "NoName_4"
    {
      Reference = "NoName_4"
      RefineWindow = rectangle [( -0.58125 1.9258 ) , ( 0.58594 2.0477 )]
    }

  # Profiles 
    Constant "ndop"
    {
      Reference = "ndop"
      EvaluateWindow
      {
        Element = rectangle [( -2 0 ) , ( 2 2 )]
        DecayLength = 0
      }
    }

    Constant "pdop"
    {
      Reference = "pdop"
      EvaluateWindow
      {
        Element = rectangle [( -0.75 1 ) , ( 0.75 2 )]
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
