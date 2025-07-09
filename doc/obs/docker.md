# just docker things

### set shm size for mpi

```shell
docker run --shm-size=3256m -w /io  --rm -it -v $PWD:/io --name phare registry.lpp.polytechnique.fr/phare/teamcity-fedora_dep:40
```
