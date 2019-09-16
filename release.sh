#! /bin/bash

state='initial'
CURR_GIT_TAGS=`git tag -l`
CURR_GIT_BRANCH=`git branch | grep \* | grep release | sed s/\*\ //`
GIT_DIRTY=`git status | egrep 'Changes to be committed|Changes not staged for commit'`
RELEASE_BRANCH=''
RELEASE_VERSION=''
while true; do
    case "$state" in
	initial)
	    if [ "$CURR_GIT_BRANCH" == "" ]; then		
		echo "Error: You are not on a release branch (e.g. release-x.y.z).  Create and FINALIZE the release branch before running this script."
		exit
	    elif [ "$GIT_DIRTY" != "" ]; then
		echo "Error: Your current branch is dirty.  Commit your changes before running this script";
		exit
	    else
		state='confirm_version'
		RELEASE_VERSION=`echo $CURR_GIT_BRANCH | sed s/release-//`

		cat <<EOF
Preparing for git version release.  Be sure your scripts are tested and documentation is up to date before running this script.

Here are the current git tags:
$CURR_GIT_TAGS

Here is your current branch: $CURR_GIT_BRANCH

It appears this release number $RELEASE_VERSION
EOF
		echo -en "Is this correct? [y/N]: ";
		state='verify'
		N_STATE='version_error'
		Y_STATE='bump_release'
	    fi
	    ;;
	
	bump_release)
	    RELEASE_BRANCH=$CURR_GIT_BRANCH
	    echo "Confirmed version $RELEASE_VERSION";
	    echo ""
	    echo "Next, ensure your ChangeLog is up to date and has the new version updates in it"
	    echo "Press [ENTER] to edit ChangeLog"
	    read -r TRASH
	    vi ChangeLog

	    echo "Ensure your manual is up to date and has the new version updates in it"
	    echo "Press [ENTER] to edit doc/manual.md"
	    read -r TRASH
	    vi doc/manual.md

	    state='verify'
	    N_STATE='bump_release'
	    Y_STATE='finalize_release'
	    echo -en "Documentation updated.\nProceed with release? [y/N]: "
	    ;;

	verify)
	    read -r VERIFY
	    if [ "$VERIFY" == 'Y' ] || [ "$VERIFY" == 'y' ]; then
		state=$Y_STATE
	    else
		state=$N_STATE
	    fi
	    ;;

	finalize_release)
	    echo -ne "Committing version bump changes\nPress [ENTER] to continue"
	    read -r TRASH
	    git commit -a -m "Version bump" || { echo "Error: Command failed" && exit 1; }

	    echo -ne "Merging develop branch\nPress [ENTER] to continue"
	    read -r TRASH
	    git checkout develop || { echo "Error: Command failed" && exit 1; }
	    git merge --no-ff $RELEASE_BRANCH || { echo "Error: Command failed" && exit 1; }

	    echo -ne "Cleaning release branch and merging with master\nPress [ENTER] to continue"
	    read -r TRASH

	    git checkout $RELEASE_BRANCH || { echo "Error: Command failed" && exit 1; }
	    git rm branch.info || { echo "Error: Command failed" && exit 1; }
	    git rm Release.checklist || { echo "Error: Command failed" && exit 1; }
	    git rm TODO || { echo "Error: Command failed" && exit 1; }
	    git rm Makefile || { echo "Error: Command failed" && exit 1; }
	    git rm release.sh || { echo "Error: Command failed" && exit 1; }
	    git commit -a -m "File removal in preparation for release" || { echo "Error: Command failed" && exit 1; }
	    git checkout master || { echo "Error: Command failed" && exit 1; }
	    git merge --no-ff $RELEASE_BRANCH || { echo "Error: Command failed" && exit 1; }

	    echo -ne "Tagging release and pushing to remote\nPress [ENTER] to continue"
	    read -r TRASH

	    git push origin || { echo "Error: Command failed" && exit 1; }
	    git tag $RELEASE_VERSION || { echo "Error: Command failed" && exit 1; }
	    git push origin --tags || { echo "Error: Command failed" && exit 1; }
	    
	    echo ""
	    echo "All operations complete.  You may delete the release branch if you desire"
	    exit
	    ;;

	version_error)
	    echo "Could not auto-detect version number from git branch name.  Fix script and run again!";
	    exit
	    ;;
	
	*)
	    echo "UNKNOWN STATE: $state";
	    exit
	    ;;
    esac
    
done
