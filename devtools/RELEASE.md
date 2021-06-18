# How to release


## Setup

Install required packages:
```
$ conda env create -f binder/environment.yml
```

Or if needed, update your conda environment:
```
$ conda env update -f binder/environment.yml
```

For Zenodo integration, see [Making Your Code Citable](https://guides.github.com/activities/citable-code/).

To publish a package in [PyPI](https://pypi.org/):

- Create an [account](https://pypi.org/account/register/).
- Create a new API token in the [account settings](https://pypi.org/manage/account/#api-tokens). Copy this token because you won't be able to see it again.
- Paste this token in the [GitHub secrets for your repo](https://docs.github.com/en/actions/reference/encrypted-secrets#creating-encrypted-secrets-for-a-repository) with the name `PYPI_API_TOKEN`.

## Tests

Before any release, double-check all tests had run successfully:
```
$ make tests
```


## Update version number

We use `bump2version` to update and synchronize the version number across different files.

For patch update (x.y.z → x.y.**z+1**):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg patch
```

For minor update (x.y.z → x.**y+1**.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg minor
```

For major update (x.y.z → **x+1**.0.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg major
```

Remark:

1. For a dry run with `bump2version`, use option `-n`.
2. `bump2version` will fail if the git working directory is not clean, i.e. all changes are not commited.

Once version number is updated, push everything to GitHub:
```
$ git push origin
$ git push origin --tags
```


## Add new release on GitHub

On [GitHub release page](https://github.com/patrickfuchs/buildH/releases) :

- Click the *Draft a release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy recent changes from `CHANGELOG.md` and paste them in the *Describe this release* field.
- Hit the *Publish Release* button :rocket:


## Zenodo integration

After the creation of the new release on GitHub, check a new archive has been created on [Zenodo](https://doi.org/10.5281/zenodo.4676217).


## PyPI package

After the creation of the new release on GitHub, check a new package has been published on [PyPI](https://pypi.org/project/buildh/).


If you need to manually build and upload your package on PyPI, run the following commands:

```bash
$ make build
$ make upload-to-pypi
```

Enter your username and password upon request.


## Bioconda package

The `meta.yaml` file used to create the conda package is under the folder `bioconda`.

Here a link to the merged PR which add buildH to the bioconda channel : https://github.com/bioconda/bioconda-recipes/pull/28673.

Once a new release is made on GitHub, it should be automatically updated on the bioconda channel (https://bioconda.github.io/contributor/updating.html).
