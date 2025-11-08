# $Id: Makefile 2031 2010-07-21 14:32:19Z maufder $

topdir := ../../../..


# we have to find the necessary include files of the module
MODULE_INCLUDES = -I..

include $(topdir)/Make.module


$(foreach lib, $(srcfiles), $(eval $(call MODULE, $(lib))))


include $(topdir)/Make.rules
