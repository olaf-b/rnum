require 'mkmf'
dir_config('lapack')
dir_config('cblas')
dir_config('blas')
dir_config('gfortran')
dir_config('g2c')

$distcleanfiles << 'ratlas_config.h'
ratlas_config = File::open('ratlas_config.h', 'w')
ratlas_config << "#ifndef RATLAS_CONFIG_H\n#define RATLAS_CONFIG_H\n"

# Checking GCC and GCC version
gcc_version = false
begin
  cc_version = `#{CONFIG['CC']} --version`
  if cc_version.scan(/\(GCC\)/).empty?
    puts("#{CONFIG['CC']} is not GCC")
  else
    gcc_version = `#{CONFIG['CC']} -dumpversion`.strip
    $CPPFLAGS += " -Wall"
  end
rescue
  puts("#{CONFIG['CC']} is not GCC")
end

with_veclib = with_config('veclib')
with_accelerate = with_config('accelerate')

cblas_h = with_config('cblas_h') || "cblas.h"
cblaslib = with_config('cblaslib') 
blaslib = with_config('blaslib')
lapacklib = with_config('lapacklib')

# If on Darwin/Mac OS, use vecLib framework - unless --with-veclib=no or
# --with-accelerate=no is given.
if RUBY_PLATFORM =~ /darwin/ and not (with_veclib || with_accelerate)
  if have_header('Accelerate/Accelerate.h')
    $LDFLAGS += " -framework Accelerate"
    ratlas_config << "#include <Accelerate/Accelerate.h>\n"
  elsif have_header('vecLib/vecLib.h')
    $LDFLAGS += " -framework vecLib"
    ratlas_config << "#include <vecLib/vecLib.h>\n"
  else
    puts("Was looking for a usable blas/lapack framework. Found none.  Giving up!")
    exit
  end
else
  if not have_header(cblas_h) and 
    not find_header(cblas_h, '/usr/local/include')
    puts("No #{cblas_h} found.  Try using --with-cblas-include=PATH")
    exit
  end
  ratlas_config << "#include <#{cblas_h}>\n"

  if cblaslib
    if not have_library('cblas')
      puts("No #{cblaslib} found.  Try using --with-cblaslib=CBLASLIB")
    end
  else
    if not have_library('atlas_r') and not have_library('atlas')
      puts("No atlas library found.  Will probably resort to reference BLAS.")
    end
    if not have_library('ptcblas') and not have_library('cblas')
      puts("No cblas library found.  Try using --with-cblas-lib=PATH")
      exit
    end
  end
  
  # Link with correct fortran libs when using gcc.
  if gcc_version
    if gcc_version < "4.0"
      if not have_library('g2c')
          puts("No g2c library found.  It might be required for lapack.")
          puts("Try using --with-g2c-lib=PATH")
      end
    else
      if not have_library('gfortran') and not find_library('gfortran',
        '', "/usr/lib/gcc-#{gcc_version}", 
        "/usr/local/lib/gcc-#{gcc_version}")
        puts("No libgfortran found.  It might be required for lapack.")
        puts("You might need to create a symlink libgfortran.so to")
        puts("libgfortran#.so, where # depends on your compiler.")
        puts("Try using --with-gfortran-lib=PATH")
      end
    end
  end

  if blaslib
    if not have_library(blaslib)
      puts("No #{blaslib} found.  Try using --with-blaslib=BLASLIB")
    end
  else
    if not have_library('ptf77blas') and not have_library('f77blas') and 
      not have_library('blas')
      puts("No fortran blas library found. It might be required for lapack.")
      puts("Try using --with-blas-lib=PATH")
    end
  end

  if lapacklib
    if not have_library(lapacklib)
      puts("No #{lapacklib} found.  Try using --with-lapacklib=LAPACKLIB")
    end
  else
    if not have_library('alapack_r') and not have_library('alapack') and 
      not have_library('lapack')
      puts("No lapack library found.  Try using --with-lapack-lib=PATH")
      exit
    end
  end
end

ratlas_config << "#endif\n"
ratlas_config.close

$INSTALLFILES = [   ['ratlas.h',        '$(archdir)'],
                    ['ratlas_config.h', '$(archdir)'],
                    ['ratlas_func.h',   '$(archdir)'],
                    ['ratlas_lapack.h', '$(archdir)'],
                    ['ratlas_cblas.h',  '$(archdir)']   ]

create_makefile("ratlas")
