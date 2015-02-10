
#-------------------------------------------------------------------------------
# Name:        test
# Purpose:     
#
# Author:      Francesca
#
# Created:     25/11/2014

#------------------------------------------------------------------------------

"""
Usage:=
"""
def testf():
    """ This is a test """
    file = open('test.txt', 'w')
    file.write("this is a test file\n")

    file.close()