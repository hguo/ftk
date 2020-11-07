#/bin/sh

mkdir -p $SPACK_ROOT &&
      git clone --depth 50 https://github.com/spack/spack.git $SPACK_ROOT &&
      printf "config:\n  build_jobs: 2\n" > $SPACK_ROOT/etc/spack/config.yaml &&
      printf "packages:\n  all:\n    target: ['x86_64']\n" \
              > $SPACK_ROOT/etc/spack/packages.yaml;

spack mirror add E4S https://cache.e4s.io
wget https://oaciss.uoregon.edu/e4s/e4s.pub
spack gpg trust e4s.pub
