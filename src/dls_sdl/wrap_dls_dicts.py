from wrap_dls_connect import *

__all__ = ['handle_enum',
            'resultcode',
            'tempmdl',
            'mwrmdl',
            'rgmdl',
            'regudist',
            'reguxval',
            'dict_tempmdl',
            'dict_mwrmdl',
            'dict_rgmdl',
            'dict_regudist',
            'dict_reguxval']

#Converts DYNAMICS enum class attribute to label string for display
#inputs: 'x', enum attribute to be converted, and 'enum_class', a variable name for a DYNAMICS enum
#outputs: dict[x], the label string corresponding to the enum attribute

#----------DYNAMICS ENUMS AS PYTHON DICTIONARIES----------

def __getattr__(name):
    dynamicsLib = get_dynamics_lib()

    resultcode = dynamicsLib.ResultCode
    tempmdl = dynamicsLib.SolventTemperatureModel
    rgmdl = dynamicsLib.RG_Model
    mwrmdl = dynamicsLib.MWR_Model
    regudist = dynamicsLib.ReguDistType
    reguxval = dynamicsLib.ReguXValues

    #Solvent temperature model
    dict_tempmdl = {
        tempmdl.SOLVTEMPMDL_NONE: 'none',
        tempmdl.SOLVTEMPMDL_FIXED: 'fixed',
        tempmdl.SOLVTEMPMDL_AQUEOUS: 'aqueous',
    }

    #Weight-averaged molar mass model
    dict_mwrmdl = {
        mwrmdl.MWR_NO_MWR_MODEL: 'no model',
        mwrmdl.MWR_GLOBULAR_PROTIENS: 'globular proteins',
        mwrmdl.MWR_LINEAR_POLYMERS: 'linear polymers',
        mwrmdl.MWR_BRANCHED_POLYMERS: 'branched polymers',
        mwrmdl.MWR_STARBURST_POLYMERS: 'starburst polymers'
    }

    #Conformation model for SLS
    dict_rgmdl = {
        rgmdl.RG_HOLLOW_SPHERE: 'hollow sphere',
        rgmdl.RG_SPHERE: 'solid sphere',
        rgmdl.RG_RANDOM_COIL: 'random coil',
        rgmdl.RG_REGULAR_STAR2ARMS: 'regular star with 2 arms',
        rgmdl.RG_REGULAR_STAR3ARMS: 'regular star with 3 arms',
        rgmdl.RG_REGULAR_STAR4ARMS: 'regular star with 4 arms',
        rgmdl.RG_REGULAR_STAR5ARMS: 'regular star with 5 arms'
    }

    #Regularization data distribution
    dict_regudist = {
        'I': regudist.RDIST_INTENSITY,
        'M': regudist.RDIST_MASS,
        'N': regudist.RDIST_NUMBER
    }

    #Regularization data x-axis quantity
    dict_reguxval = {
        'C': reguxval.RXTYPE_DIFF_COEFF,
        'D': reguxval.RXTYPE_DIAMETER,
        'R': reguxval.RXTYPE_RADIUS,
        'T': reguxval.RXTYPE_DECAY_TIME
    }

    def handle_enum(x, enum_class):
        dict={}
        if enum_class==tempmdl: dict = dict_tempmdl
        elif enum_class==mwrmdl: dict = dict_mwrmdl
        elif enum_class==rgmdl: dict = dict_rgmdl
        elif enum_class==regudist: dict = dict_regudist
        elif enum_class==reguxval: dict = dict_reguxval
        return dict[x]

    mapping = {
        'resultcode': resultcode,
        'tempmdl': tempmdl,
        'mwrmdl': mwrmdl,
        'rgmdl': rgmdl,
        'regudist': regudist,
        'reguxval': reguxval,
        'dict_tempmdl': dict_tempmdl,
        'dict_mwrmdl': dict_mwrmdl,
        'dict_rgmdl': dict_rgmdl,
        'dict_regudist': dict_regudist,
        'dict_reguxval': dict_reguxval,
        'handle_enum': handle_enum
    }

    if name in mapping:
        return mapping[name]
    
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
