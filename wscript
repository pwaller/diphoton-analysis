#! /usr/bin/env python

from waflib.TaskGen import feature, after
from waflib.Task import Task
from waflib.Tools import c_preproc

@feature('cxx')
@after('apply_link')
def process_pch(self):
	if getattr(self, 'pch', ''):
		nodes = self.to_nodes(self.pch)
		for x in nodes:
			self.create_task('gchx', x, x.change_ext('.h.gch'))

class gchx(Task):
	run_str = ('${CXX} ${CXXFLAGS} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} '
	           '${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F}${SRC}' 
	           ' ${CXX_TGT_F}${TGT}')
	scan    = c_preproc.scan
	ext_out = ['.h']
	color   = 'BLUE'

from waflib.Task import update_outputs
update_outputs(gchx)

def options(opt):
    opt.load('compiler_c compiler_cxx python')
    opt.load('proto check_with', tooldir="./common/waf")
    
    opt.add_option('--with-a4', default=None,
        help="Also look for a4 at the given path")

def configure(conf):
    conf.load('compiler_c compiler_cxx python')
    conf.load('proto check_with', tooldir="./common/waf")
    
    conf.check_with(conf.check_cfg, "a4", package="a4", args="--cflags --libs")
    
    #conf.check_cfg(package="a4", uselib_store="A4", args="--libs --cflags")
    conf.env.RPATH_A4 = conf.env.LIBPATH_A4
    conf.env.LIB_A4.insert(0, "a4atlas")
    #conf.env.STLIB_A4.insert(0, "a4atlas")
    #conf.env.STLIB_A4.insert(0, "c")
    #conf.env.STLIB_A4.append("rt")
    #conf.env.STLIB_A4.append("dl")
    #conf.env.STLIB_A4.append("c")
    
    conf.env.append_value("CXXFLAGS", ["-std=c++0x", "-ggdb"])
    conf.env.append_value("LDFLAGS", ["-Wl,--as-needed"])
    conf.env.append_value("RPATH", [conf.env.LIBDIR])
    
    conf.check_cfg(path="root-config", package="", uselib_store="CERN_ROOT_SYSTEM",
                   args='--libs --cflags', mandatory=False)
                   
    conf.to_log("Final environment:")
    conf.to_log(conf.env)

def build(bld):
    bld.shlib(features="cxx cxxstlib", source=bld.path.ant_glob("proto/**/*.proto"), 
              target="analysis_protobuf", use=["A4", ""])

    bld.shlib(features="cxxstlib", source=bld.path.ant_glob("src/external/**.cxx"), 
              target="analysis_externals", use=["CERN_ROOT_SYSTEM"])
    
    bld.program(
        source=bld.path.ant_glob("src/**.cxx"),
        target="analysis",
        includes="pch . src src/external",
        pch="pch/all.h",
        use=["analysis_externals", "analysis_protobuf", "A4"],
    )
    
    for path in bld.path.ant_glob("src/apps/**.cxx"):
        bld.program(
            "cxx",
            source=[path],
            includes="pch src src/external",
            target=path.name[:-len(".cxx")],
            use=["analysis_externals", "analysis_protobuf", "A4"],
        )
