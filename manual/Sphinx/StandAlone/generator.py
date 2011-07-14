import os
import string

# variables

#----Folders

gettingstarted = "GettingStarted"
theory   = "Theory"
referenceguide = "ReferenceGuide"
examples = "Examples"

#----filenames

inputfile = "InputFile.rst" # only for Getting Started

elasticity = "Elasticity.rst"
driftdiffusion = "DriftDiffusion.rst"
thermal    = "Thermal.rst"
efa        = "EFA.rst"
dsc        = "Dsc.rst"


#----Output Folders 

freedownload   = "FreeDownload"
standalone =     "StandAlone"
multiusers =     "MultiUsers"


# function that copies data between any couple of markers


def checkMarker(toControl, toExport) :
     
     if not(os.path.isfile(toControl)) :
        return "\n\t--->" + toControl + " is not present\n" # is the file present?
     
     checking = open(toControl, "r")
     copied = open( toExport, "a")
     
     counter = 0 # for <marker>
     ending = 0  # for </marker>
     row = checking.readline()
     while len(row) > 0 :

              if len(row) == 0 :
                break
              if row.find("</marker>") > -1 :
                ending = ending + 1
              if row.find("<marker>") > -1 :
                    counter = counter + 1
                    while row.find("</marker>") == -1 and len(row) > 0 :
                        row = checking.readline()
                        if row.find("</marker>") > -1 :
                            break
                        copied.write(row)
                        
                    ending = ending + 1
              row = checking.readline()

     if counter != ending :
        return "\n\nATTENTION :    Inequality on markers detected!" + "\n\ton file : " + toControl + "\nPlease, check this file!"
     
     checking.close()
     copied.close()
     return "\nI've checked " + str(toControl) + ". \n\n There's parity of " + str(counter) + " couple of markers"

# where am I

currentDir = os.getcwd() #returning working directory
Folder = os.path.split(currentDir)
Foldername = Folder[1]  # only the name of current directory
     
     
# preparing paths
sourceGet = os.path.join(Folder[0], gettingstarted)
sourceTheo = os.path.join(Folder[0], theory)
sourceRef = os.path.join(Folder[0], referenceguide)
    
# Reset files

def  resetFiles(toReset) :
    if os.path.isfile(os.path.join(currentDir, toReset)) :
        os.remove(os.path.join(currentDir, toReset))
    else :
        return toReset + " is not present"

    return toReset + "prepared.\n"
    
# Let's roll

if Foldername == freedownload :
    print resetFiles("GettingStarted.rst")
    
    #  GettingStarted folder
    print checkMarker(os.path.join(sourceGet, inputfile), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, elasticity), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, driftdiffusion), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, thermal), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, efa), "GettingStarted.rst")  
    print checkMarker(os.path.join(sourceGet, dsc), "GettingStarted.rst")

    print "I've made the FreeDownload book" # make the FreeDownload book


elif Foldername == standalone :
    print resetFiles("GettingStarted.rst")
    print resetFiles("Theory.rst")
    
    #  GettingStarted folder
    print checkMarker(os.path.join(sourceGet, inputfile), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, elasticity), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, driftdiffusion), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, thermal), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, efa), "GettingStarted.rst")  
    print checkMarker(os.path.join(sourceGet, dsc), "GettingStarted.rst")
    
    #  Theory folder
    print checkMarker(os.path.join(sourceTheo, elasticity), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, driftdiffusion), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, thermal), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, efa), "Theory.rst")  
    print checkMarker(os.path.join(sourceTheo, dsc), "Theory.rst")
    
    print "I've made the StandAlone book" # make the StandAlone book


elif Foldername == multiusers :
    print resetFiles("GettingStarted.rst")
    print resetFiles("Theory.rst")
    print resetFiles("Examples.rst")
    
    #  GettingStarted folder
    print checkMarker(os.path.join(sourceGet, inputfile), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, elasticity), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, driftdiffusion), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, thermal), "GettingStarted.rst")
    print checkMarker(os.path.join(sourceGet, efa), "GettingStarted.rst")  
    print checkMarker(os.path.join(sourceGet, dsc), "GettingStarted.rst")
    
    #  Theory folder
    print checkMarker(os.path.join(sourceTheo, elasticity), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, driftdiffusion), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, thermal), "Theory.rst")
    print checkMarker(os.path.join(sourceTheo, efa), "Theory.rst")  
    print checkMarker(os.path.join(sourceTheo, dsc), "Theory.rst")
    
    #  Examples folder
    print checkMarker(os.path.join(sourceRef, elasticity), "Examples.rst")
    print checkMarker(os.path.join(sourceRef, driftdiffusion), "Examples.rst")
    print checkMarker(os.path.join(sourceRef, thermal), "Examples.rst")
    print checkMarker(os.path.join(sourceRef, efa), "Examples.rst")  
    print checkMarker(os.path.join(sourceRef, dsc), "Examples.rst")

    print "I've made the MultiUsers book" # make the MultiUsers book
    
else :
    print "You are in a wrong directory!"

# Oh No!, I've just finished


