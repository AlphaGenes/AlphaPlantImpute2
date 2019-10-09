import setuptools

def our_version():
    """Customise the version string from git tags and commits.
    By default setuptools_scm gives a highly detailed version string. We want a
    simplified one that is just the version or version plus the commit if
    there are newer commits. For example:
      Tag is the most recent commit      -  0.1.2
      Commits after the most recent tag  -  0.1.2+gbccb995
    """
    def version_scheme(version):
        if version.distance is None or version.distance == 0:
            # No commits after most recent tag - display tag only
            return version.format_with("{tag}")
        else:
            # Commits after most recent tag - display tag and most recent commit
            return version.format_with("{tag}+{node}")

    def local_scheme(version):
        # Ignore dirty/clean status of git repo completely
        return ''

    return {'version_scheme': version_scheme, 'local_scheme': local_scheme}


setuptools.setup(
    use_scm_version=our_version
)
