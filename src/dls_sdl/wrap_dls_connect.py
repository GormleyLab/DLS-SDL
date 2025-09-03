import sys
import clr

dynamics_path = r'C:\Program Files (x86)\WTC\DYNAMICS SDK 3.0.1\Src\ComTest\Bin\Release'

dynamicsLib = None
dynapro = None

#----------CONNECT------------
def check_deps():
    test = Assembly.GetReferencedAssemblies(clr.AddReference(r'DynamicsLib'))
    print('Dependencies:')
    [print(i) for i in test]
    print('\nLoaded:')
    [print(i) for i in clr.ListAssemblies(False)]

def connect_to_dls_api(path):
    clr.AddReference("System")
    #from System import UInt32, Int32
    from System.Reflection import Assembly

    #Add dynamicsLib.dll assembly file to path and check dependencies
    sys.path.append(path)
    #check_deps()

    #Connect to the dynamicsLib.dll assembly file
    clr.AddReference(r'DynamicsLib')
    import dynamicsLib
    return dynamicsLib

def connect_to_dls_instrument(dynamicsLib):
    return dynamicsLib.CControlServer()

def get_dynamics_lib():
    global dynamicsLib
    global dynamics_path
    if dynamicsLib is None:
        dynamicsLib = connect_to_dls_api(dynamics_path)
    return dynamicsLib

def get_dynapro(dynamicsLib):
    global dynapro
    if dynapro is None:
        dynapro = connect_to_dls_instrument(dynamicsLib)
    return dynapro