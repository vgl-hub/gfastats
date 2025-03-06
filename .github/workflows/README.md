validate is automatically run on pushes to any branch, or pull requests to main

to automatically create a new release and automatically upload mac, ubuntu, and windows builds run:
`git tag v*`
`git push origin v*`
where * is the version number.

Example:
`git tag v1.2.1`
`git push origin v1.2.1`

