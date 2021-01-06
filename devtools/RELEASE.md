# How to release


## Setup

Install required packages:
```
$ conda create env -f binder/environment.yml
```

If needed, update your conda environment:
```
$ conda env update -f binder/environment.yml
```

## Tests

Before any release, double-check all tests had run successfully:
```
$ make tests
```


## Update version number

We use `bump2version` to update and synchronize the version number across different files.

For patch update (x.y.**z**):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg patch
```

For minor update (x.**y**.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg minor
```

For major update (**x**.0.0):
```
$ bump2version --verbose --config-file devtools/bumpversion.cfg major
```

Remarque : 

1. For a dry run with `bump2version`, use option `-n`.
2. `bump2version` will fail if the git working directory is not clean, i.e. all changes are not commited.

Once version number is updated, push everything to GitHub:
```
$ git push origin
$ git push origin --tags
```


### Add new release on GitHub

On [GitHub release page](https://github.com/pierrepo/autoclasswrapper/releases) :

- Click the *Draft a release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy recent changes from `CHANGELOG.md` and paste them in the *Describe this release* field.
- Hit the *Publish Release* button :rocket:

