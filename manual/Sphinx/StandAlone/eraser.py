import os
import shutil
import string

# imposing variables

gettingstarted = "GettingStarted.rst"
theory   = "Theory.rst"
referenceguide = "ReferenceGuide.rst"
examples =       "Examples.rst"
manual   = "TiberCAD.pdf" # it's important that this name is right!

# imposing Folder

freedownload   = "FreeDownload"
standalone =     "StandAlone"
multiusers =     "MultiUsers"

currentDir = os.getcwd() #returning working directory

Folder = os.path.split(currentDir)
Foldername = Folder[1]  # only the name of current directory

# Move result in Sphinx folder

buildpath = os.path.join(currentDir, "_build")
latexpath = os.path.join(buildpath, "latex" )

if Foldername == freedownload :
    shutil.copy(os.path.join(latexpath, manual), os.path.join(currentDir, "TiberCAD-FreeDownload.pdf")) # thanks shell utils!

elif Foldername == standalone :
    shutil.copy(os.path.join(latexpath, manual), os.path.join(currentDir, "TiberCAD-GettingStarted.pdf")) # thanks shell utils!

elif Foldername == multiusers :
    shutil.copy(os.path.join(latexpath, manual), os.path.join(currentDir, "TiberCAD-MultiUsers.pdf")) # thanks shell utils!

else :
    print "You are in a wrong directory. \n\nBe careful you've just attempted to delete this folder!"
    

# I'll erase every .rst and the _build folder here

# delete .rst

if os.path.isfile(os.path.join(currentDir, gettingstarted)) :
    os.remove(os.path.join(currentDir, gettingstarted))
    print "\n\t" + str(gettingstarted) + " removed."

if os.path.isfile(os.path.join(currentDir, theory)) :
    os.remove(os.path.join(currentDir, theory))
    print "\n\t" + str(theory) + " removed."
    
if os.path.isfile(os.path.join(currentDir, referenceguide)) :
    os.remove(os.path.join(currentDir, referenceguide))
    print "\n\t" + str(referenceguide) + " removed."
    
if os.path.isfile(os.path.join(currentDir, examples)) :
    os.remove(os.path.join(currentDir, examples))
    print "\n\t" + str(examples) + " removed."


# delete temp folder -BEWARE

if os.path.isdir(os.path.join(currentDir, "_build")) and ( Foldername == "MultiUsers" or Foldername == "StandAlone" or Foldername == "FreeDownload") :
    shutil.rmtree(os.path.join(currentDir, "_build"))
    print "\n\n\t--->_Build Folder removed."
