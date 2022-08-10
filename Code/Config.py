import os,pickle,sys

if os.uname()[1].split('.')[1]=='frontera':
    config = {
    'storm':{
    'z0':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.004096'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096.z0.000.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096.z0.000.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096.z0.000.AHF_halos',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536.z0.000.AHF_halos',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.004096.z0.000.AHF_halos'
        }
    },
    'z0.04':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.003936',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.003936',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.003936',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.003936',
            'old_vdXsec':'',
            'vdXsec':''
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.003936.z0.041.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.003936.z0.041.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.003936.z0.041.AHF_halos',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.003936.z0.041.AHF_halos',
            'old_vdXsec':'',
            'vdXsec':''
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.001813',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.001813',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.001813',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.028992',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.001813'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.001813.z0.999.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.001813.z0.999.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.001813.z0.999.AHF_halos',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.028992.z1.000.AHF_halos',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.001813.z0.999.AHF_halos'
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.001025',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.001025',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.001025',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.001056',
            'old_vdXsec':'',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.001025'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.001025.z1.999.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.001025.z1.999.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.001025.z1.999.AHF_halos',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.001056.z1.938.AHF_halos',
            'old_vdXsec':'',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.001025.z1.999.AHF_halos'
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000672',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000672',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000672',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000672',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.010560',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000672'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000672.z2.998.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000672.z2.998.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000672.z2.998.AHF_halos',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000672.z2.998.AHF_halos',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.010560.z3.047.AHF_halos',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000672.z2.998.AHF_halos'
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000482',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000482',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000482',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000480',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.007703',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000482.z3'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000482.z3.996.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000482.z3.996.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000482.z3.996.AHF_halos',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000480.z4.010.AHF_halos',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.007703.z4.000.AHF_halos',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000482.z3.996.AHF_halos'
        }
    },
    'z10':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000146',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000146',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000146',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.002334',
            'vdXsec':''
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000146.z9.991.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000146.z9.991.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000146.z9.991.AHF_halos',
            'SI50':'',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.002334.z9.997.AHF_halos',
            'vdXsec':''
        }
    },
    'z13.5':{
        'simpaths':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000096',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000096',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000096',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000096',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.001536',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000096'
        },
        'AHFs':{
            'CDM':'/scratch1/08902/tg882017/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.000096.z13.428.AHF_halos',
            'SI3':'/scratch1/08902/tg882017/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.000096.z13.428.AHF_halos',
            'SI10':'/scratch1/08902/tg882017/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.000096.z13.428.AHF_halos',
            'SI50':'/scratch1/08902/tg882017/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.000096.z13.428.AHF_halos',
            'old_vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.001536.z13.428.AHF_halos',
            'vdXsec':'/scratch1/08902/tg882017/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.000096.z13.428.AHF_halos'
        }
    }
    },


    'h148':{
    'z0':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':''
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.027836'
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.027836.z1.000.AHF_halos'
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.015580'
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.015580.z2.000.AHF_halos'
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.010182'
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'/scratch1/08902/tg882017/h148.cosmo50PLKvdXsec.3072/h148.cosmo50PLKvdXsec.3072.010182.z3.000.AHF_halos'
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'',
            'h148':''
        },
        'AHFs':{
            'CDM':'',
            'h148':''
        }
    }
    },


    'gitdir':'/work2/08902/tg882017/SelfInteractingDarkMatter',
    'python':'/home1/08902/tg882017/anaconda3/bin/python'
    }




elif os.uname()[1].split('.')[0]=='glaurung':
    config = {
    'storm':{
    'z0':{
        'simpaths':{
            'CDM':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096',
            'SI3':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096',
            'SI10':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096',
            'SI30':'',
            'SI50':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.003936',
            'old_vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536',
            'vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.004096'
        },
        'AHFs':{
            'CDM':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmb.4096/storm.cosmo25cmb.4096.004096.z0.000.AHF_halos',
            'SI3':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI3.4096/storm.cosmo25cmbSI3.4096.004096.z0.000.AHF_halos',
            'SI10':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI10.4096/storm.cosmo25cmbSI10.4096.004096.z0.000.AHF_halos',
            'SI30':'',
            'SI50':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbSI50.4096/storm.cosmo25cmbSI50.4096.003936.z0.041.AHF_halos',
            'old_vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.065536.z0.000.AHF_halos',
            'vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.4096.VTS/storm.cosmo25cmbvdXsec.4096.VTS.004096.z0.000.AHF_halos'
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.028992',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'/home/vannest/dwarf_volumes/storm.SIDM/storm.cosmo25cmbvdXsec.65536/storm.cosmo25cmbvdXsec.65536.028992.z1.000.AHF_halos',
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    }
    },


    'h148':{
    'z0':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'',
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'',
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'',
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'',
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'old_vdXsec':'',
        }
    }
    },


    'gitdir':'/myhome2/users/vannest/SelfInteractingDarkMatter',
    'python':'/myhome2/users/vannest/anaconda3/bin/python'
    }





elif os.uname()[1].split('.')[0]=='dodo':
    config = {
    'storm':{
    'z0':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'',
            'SI3':'',
            'SI10':'',
            'old_vdXsec':'',
        }
    }
    },


    'h148':{
    'z0':{
        'simpaths':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.004096',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.004096.0000.z0.000.AHF_halos',
            'old_vdXsec':'',
        }
    },
    'z1':{
        'simpaths':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.001740',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.001740.0000.z1.000.AHF_halos',
            'old_vdXsec':'',
        }
    },
    'z2':{
        'simpaths':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000974',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000974.0000.z1.999.AHF_halos',
            'old_vdXsec':'',
        }
    },
    'z3':{
        'simpaths':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000633',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000633.0000.z3.014.AHF_halos',
            'old_vdXsec':'',
        }
    },
    'z4':{
        'simpaths':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000475',
            'old_vdXsec':'',
        },
        'AHFs':{
            'CDM':'/data/REPOSITORY/e12Gals/h148.cosmo50PLK.3072/h148.cosmo50PLK.3072.000475.0000.z3.864.AHF_halos',
            'old_vdXsec':'',
        }
    }
    },


    'gitdir':'/home/jv590/SelfInteractingDarkMatter',
    'python':'/home/jv590/anaconda3/bin/python'
    }



else:
    print('Check config settings')
    sys.exit(0)

out = open('Config.pickle','wb')
pickle.dump(config,out)
out.close()
