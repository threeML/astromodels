cd BUILD_DIR

# We need a custom include and library path to use the packages installed
# in the build environment
export CFLAGS="${CFLAGS} -I${PREFIX}/include -Wall --pedantic -Wno-comment -Wno-long-long -g  -ffloat-store -fPIC"
export CXXFLAGS="${CXXFLAGS} -I${PREFIX}/include -Wall --pedantic -Wno-comment -Wno-long-long -g  -ffloat-store -fPIC -std=c++14"
export CPPFLAGS="${CXXFLAGS} -I${PREFIX}/include"
export LDFLAGS="${LDFLAGS} -L${PREFIX}/lib"

if [ "$(uname)" == "Linux" ]; then
    
    ./configure --prefix=${SRC_DIR}/xspec-modelsonly-build
    
    ./hmake 'XSLM_USER_FLAGS="-I${PREFIX}/include"' 'XSLM_USER_LIBS="-L${PREFIX}/lib -lCCfits -lcfitsio -lwcslib -lgfortran"'

fi

if [ "$(uname)" == "Darwin" ]; then
    
    # Build for a fairly old mac to ensure portability
        
    ./configure --prefix=${SRC_DIR}/xspec-modelsonly-build

    ./hmake 'LDFLAGS_CXX=-headerpad_max_install_names -lcfitsio -lCCfits -lccfits -lwcs -lgfortran' 'XSLM_USER_LIBS="-L${PREFIX}/lib -lCCfits -lcfitsio -lwcslib -lgfortran"'

fi

make install

# Correct the output of make install so that the output directory structure
# makes more sense for a conda environment. We will place libraries in the
# ${PREFIX}/lib directory and data files in ${PREFIX}/Xspec/spectral

# Copy libraries in the ${PREFIX}/lib directory

cp -v `find ${SRC_DIR}/xspec-modelsonly-build/Xspec/ -name "libXS*"` ${PREFIX}/lib

# Create a ${PREFIX}/lib/Xspec directory
mkdir ${PREFIX}/lib/Xspec

# Create an empty headas directory. The env. variable HEADAS should point here at
# runtime
mkdir ${PREFIX}/lib/Xspec/headas

# Fill it with a useless file, otherwise Conda will remove it during installation
echo "LEAVE IT HERE" > ${PREFIX}/lib/Xspec/headas/DO_NOT_REMOVE

# Copy the spectral data in the right place. According to the Xspec documentation,
# this should be $HEADAS/../spectral
cp -rv ${SRC_DIR}/xspec-modelsonly-build/spectral ${PREFIX}/lib/Xspec
