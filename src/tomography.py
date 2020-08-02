# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_tomography')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_tomography')
    _tomography = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_tomography', [dirname(__file__)])
        except ImportError:
            import _tomography
            return _tomography
        try:
            _mod = imp.load_module('_tomography', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _tomography = swig_import_helper()
    del swig_import_helper
else:
    import _tomography
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class dp(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, dp, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, dp, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _tomography.new_dp()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _tomography.delete_dp
    __del__ = lambda self: None

    def assign(self, value):
        return _tomography.dp_assign(self, value)

    def value(self):
        return _tomography.dp_value(self)

    def cast(self):
        return _tomography.dp_cast(self)
    if _newclass:
        frompointer = staticmethod(_tomography.dp_frompointer)
    else:
        frompointer = _tomography.dp_frompointer
dp_swigregister = _tomography.dp_swigregister
dp_swigregister(dp)

def dp_frompointer(t):
    return _tomography.dp_frompointer(t)
dp_frompointer = _tomography.dp_frompointer

class tomography(_object):
    """Proxy of C++ tomography class."""

    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, tomography, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, tomography, name)
    __repr__ = _swig_repr
    __swig_setmethods__["tree"] = _tomography.tomography_tree_set
    __swig_getmethods__["tree"] = _tomography.tomography_tree_get
    if _newclass:
        tree = _swig_property(_tomography.tomography_tree_get, _tomography.tomography_tree_set)

    def __init__(self, *args):
        """
        __init__(tomography self) -> tomography
        __init__(tomography self, double xrmin, double xrmax, int xrbins, double xzsrcmin, double xzsrcmax, char * xoutfile, bool xverbose) -> tomography
        """
        this = _tomography.new_tomography(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _tomography.delete_tomography
    __del__ = lambda self: None

    def allocate_lens_memory(self, xNcen):
        """allocate_lens_memory(tomography self, int xNcen) -> int"""
        return _tomography.tomography_allocate_lens_memory(self, xNcen)


    def process_lens(self, xra, xdec, xzred, xwt):
        """process_lens(tomography self, double xra, double xdec, double xzred, double xwt) -> int"""
        return _tomography.tomography_process_lens(self, xra, xdec, xzred, xwt)


    def finalize_lenses(self):
        """finalize_lenses(tomography self) -> int"""
        return _tomography.tomography_finalize_lenses(self)


    def process_source(self, sra, sdec, se1, se2, swt, smcat, sc1_dp, sc2_dp, sc1_nb, sc2_nb, szbest, usepdf):
        """process_source(tomography self, double sra, double sdec, double se1, double se2, double swt, double smcat, double sc1_dp, double sc2_dp, double sc1_nb, double sc2_nb, double szbest, bool usepdf) -> int"""
        return _tomography.tomography_process_source(self, sra, sdec, se1, se2, swt, smcat, sc1_dp, sc2_dp, sc1_nb, sc2_nb, szbest, usepdf)


    def setup_pofz(self, pofz_zmin, pofz_zdiff, xNpz):
        """setup_pofz(tomography self, double pofz_zmin, double pofz_zdiff, int xNpz) -> int"""
        return _tomography.tomography_setup_pofz(self, pofz_zmin, pofz_zdiff, xNpz)


    def process_pofz(self, pofz):
        """process_pofz(tomography self, double * pofz) -> int"""
        return _tomography.tomography_process_pofz(self, pofz)


    def finalize_results(self):
        """finalize_results(tomography self) -> int"""
        return _tomography.tomography_finalize_results(self)


    def test_searchrecord(self):
        """test_searchrecord(tomography self) -> int"""
        return _tomography.tomography_test_searchrecord(self)

tomography_swigregister = _tomography.tomography_swigregister
tomography_swigregister(tomography)

# This file is compatible with both classic and new-style classes.


