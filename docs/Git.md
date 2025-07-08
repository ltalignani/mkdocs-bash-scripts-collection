# Git Version Control Guide

## Introduction

**Git** is a distributed version control system that tracks changes in files and coordinates work among multiple people. Think of it as a sophisticated "save" system that keeps track of every change you make to your project files, allowing you to go back to previous versions, compare changes, and collaborate with others safely.

This guide will teach you Git from the ground up, assuming you've never used version control before.

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Why Use Git?

Imagine working on a project and accidentally deleting important code, or wanting to try a risky change without losing your working version. Git solves these problems by:

- **Tracking all changes** to your files over time
- **Allowing you to revert** to any previous version
- **Enabling collaboration** with multiple developers
- **Backing up your work** across multiple locations
- **Managing different versions** of your project simultaneously

## Key Concepts

### Repository (Repo)
A repository is a directory that contains your project files and the complete history of changes. Think of it as a special folder that Git monitors.

### Commit
A commit is a snapshot of your project at a specific point in time. Each commit has a unique identifier and contains information about what changed.

### Branch
A branch is a parallel version of your project. You can work on different features in separate branches without affecting the main code.

### Working Directory, Staging Area, and Repository

Git has three main areas where your files can exist:

```
Working Directory → Staging Area → Repository
     (edit)           (add)        (commit)
```

- **Working Directory**: Where you edit your files
- **Staging Area**: Where you prepare changes before committing
- **Repository**: Where Git stores committed changes permanently

## Installation and Setup

### Installing Git

**On Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install git
```

**On macOS:**
```bash
# Using Homebrew
brew install git

# Or download from https://git-scm.com/
```

**On Windows:**
Download from [https://git-scm.com/download/win](https://git-scm.com/download/win)

### Initial Configuration

Before using Git, set up your identity:

```bash
# Set your name (replace with your actual name)
git config --global user.name "Your Name"

# Set your email (replace with your actual email)
git config --global user.email "your.email@example.com"

# Set default text editor (optional)
git config --global core.editor "nano"

# Check your configuration
git config --list
```

## Basic Git Workflow

### 1. Creating a Repository

#### Option A: Start a new project
```bash
# Create a new directory
mkdir my-project
cd my-project

# Initialize Git repository
git init
```

#### Option B: Clone an existing project
```bash
# Clone from a remote repository
git clone https://github.com/username/repository-name.git
cd repository-name
```

### 2. Check Repository Status

```bash
# See what's happening in your repository
git status
```

**Example output:**
```
On branch main
No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	hello.txt

nothing added to commit but untracked files present (use "git add" to track)
```

### 3. Add Files to Staging Area

```bash
# Add a specific file
git add filename.txt

# Add all files in current directory
git add .

# Add all files with specific extension
git add *.py

# Add multiple specific files
git add file1.txt file2.py file3.html
```

**Visual representation:**
```
Working Directory    Staging Area    Repository
┌─────────────────┐  ┌─────────────┐  ┌─────────────┐
│ hello.txt       │  │             │  │             │
│ (modified)      │  │             │  │             │
└─────────────────┘  └─────────────┘  └─────────────┘
                              ↓
                        git add hello.txt
                              ↓
┌─────────────────┐  ┌─────────────┐  ┌─────────────┐
│ hello.txt       │  │ hello.txt   │  │             │
│ (modified)      │  │ (staged)    │  │             │
└─────────────────┘  └─────────────┘  └─────────────┘
```

### 4. Commit Changes

```bash
# Commit with a message
git commit -m "Add hello.txt file"

# Commit with a detailed message
git commit -m "Add user authentication feature

- Added login form
- Implemented password validation
- Created user session management"

# Add and commit in one step (for tracked files only)
git commit -am "Quick commit message"
```

**After committing:**
```
Working Directory    Staging Area    Repository
┌─────────────────┐  ┌─────────────┐  ┌─────────────┐
│ hello.txt       │  │             │  │ hello.txt   │
│ (clean)         │  │             │  │ (committed) │
└─────────────────┘  └─────────────┘  └─────────────┘
```

### 5. View Commit History

```bash
# View commit history
git log

# View compact history
git log --oneline

# View history with graph
git log --oneline --graph

# View last 5 commits
git log -5
```

**Example output:**
```
commit a1b2c3d4e5f6g7h8i9j0k1l2m3n4o5p6q7r8s9t0
Author: Your Name <your.email@example.com>
Date:   Mon Oct 23 14:30:00 2023 +0200

    Add hello.txt file

commit b2c3d4e5f6g7h8i9j0k1l2m3n4o5p6q7r8s9t0a1
Author: Your Name <your.email@example.com>
Date:   Mon Oct 23 14:25:00 2023 +0200

    Initial commit
```

## Working with Files

### Checking File Status

```bash
# Detailed status
git status

# Short status
git status -s
```

**Status symbols:**  
- `??` = Untracked file  
- `A` = Added to staging area  
- `M` = Modified  
- `D` = Deleted  
- `R` = Renamed. 

### Viewing Changes

```bash
# See changes in working directory (not staged)
git diff

# See changes in staging area (staged but not committed)
git diff --staged

# See changes in a specific file
git diff filename.txt

# Compare two commits
git diff commit1 commit2
```

### Undoing Changes

```bash
# Discard changes in working directory
git checkout -- filename.txt

# Remove file from staging area (unstage)
git reset HEAD filename.txt

# Undo last commit (keep changes in working directory)
git reset --soft HEAD~1

# Undo last commit (discard changes completely)
git reset --hard HEAD~1
```

**⚠️ Warning:** `git reset --hard` permanently deletes changes!

## Branching and Merging

### Understanding Branches

Branches allow you to work on different features simultaneously:

```
main branch:     A---B---C---D
                      \
feature branch:        E---F
```

### Branch Commands

```bash
# List all branches
git branch

# Create a new branch
git branch feature-login

# Switch to a branch
git checkout feature-login

# Create and switch to a branch in one command
git checkout -b feature-login

# Delete a branch
git branch -d feature-login

# Force delete a branch
git branch -D feature-login
```

### Working with Branches

**Example workflow:**
```bash
# 1. Create and switch to feature branch
git checkout -b add-user-profile

# 2. Make changes and commit
echo "User profile page" > profile.html
git add profile.html
git commit -m "Add user profile page"

# 3. Switch back to main branch
git checkout main

# 4. Merge feature branch into main
git merge add-user-profile

# 5. Delete feature branch (optional)
git branch -d add-user-profile
```

### Merge Conflicts

When Git can't automatically merge changes, you'll get a conflict:

```bash
# This might create a conflict
git merge feature-branch
```

**Conflict markers in files:**
```
<<<<<<< HEAD
This is the content from the current branch
=======
This is the content from the feature branch
>>>>>>> feature-branch
```

**To resolve:**  
1. Edit the file to choose which content to keep  
2. Remove the conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`)  
3. Add and commit the resolved file  

```bash
# After resolving conflicts
git add conflicted-file.txt
git commit -m "Resolve merge conflict"
```

## Remote Repositories

### Understanding Remotes

A remote repository is a version of your project hosted on a server (like GitHub, GitLab, or Bitbucket).

```
Local Repository   ←→  Remote Repository
(your computer)         (GitHub/GitLab)
```

### Working with Remotes

```bash
# Add a remote repository
git remote add origin https://github.com/username/repository.git

# List remote repositories
git remote -v

# Push changes to remote
git push origin main

# Pull changes from remote
git pull origin main

# Fetch changes without merging
git fetch origin

# Clone a remote repository
git clone https://github.com/username/repository.git
```

### Typical Remote Workflow

```bash
# 1. Clone a repository
git clone https://github.com/username/project.git
cd project

# 2. Create a feature branch
git checkout -b new-feature

# 3. Make changes and commit
echo "New feature" > feature.txt
git add feature.txt
git commit -m "Add new feature"

# 4. Push branch to remote
git push origin new-feature

# 5. Switch to main and update
git checkout main
git pull origin main

# 6. Merge feature branch
git merge new-feature

# 7. Push updated main branch
git push origin main
```

## Practical Examples

### Example 1: Starting a New Project

```bash
# Create project directory
mkdir my-website
cd my-website

# Initialize Git
git init

# Create initial files
echo "# My Website" > README.md
echo "Hello, World!" > index.html

# Add and commit
git add .
git commit -m "Initial commit with README and index.html"

# Check status
git status
git log --oneline
```

### Example 2: Collaborating on a Project

```bash
# Clone existing project
git clone https://github.com/team/project.git
cd project

# Create feature branch
git checkout -b fix-bug-123

# Make changes
echo "Bug fix code" >> bugfix.py
git add bugfix.py
git commit -m "Fix bug #123: Resolve login issue"

# Push to remote
git push origin fix-bug-123

# Switch to main and update
git checkout main
git pull origin main

# Merge fix
git merge fix-bug-123
git push origin main

# Clean up
git branch -d fix-bug-123
```

### Example 3: Recovering from Mistakes

```bash
# Oops! Made a mistake in the last commit
git log --oneline
# abc123 Wrong commit message
# def456 Previous good commit

# Fix the commit message
git commit --amend -m "Correct commit message"

# Or undo the commit entirely
git reset --soft HEAD~1

# Discard all changes in working directory
git checkout -- .
```

## Best Practices

### 1. Write Good Commit Messages

**Good:**
```bash
git commit -m "Add user authentication system

- Implement login and logout functionality
- Add password encryption
- Create user session management
- Add input validation for forms"
```

**Bad:**
```bash
git commit -m "stuff"
git commit -m "fixes"
git commit -m "update"
```

### 2. Commit Often, Push Regularly

```bash
# Good practice: small, frequent commits
git commit -m "Add login form HTML"
git commit -m "Add CSS styling for login form"
git commit -m "Add JavaScript form validation"

# Then push when ready
git push origin main
```

### 3. Use Branches for Features

```bash
# Create branches for each feature
git checkout -b user-registration
git checkout -b password-reset
git checkout -b email-notifications
```

### 4. Keep Main Branch Clean

```bash
# Always test before merging to main
git checkout feature-branch
# ... test your changes ...
git checkout main
git merge feature-branch
```

## Common Commands Reference

### Essential Commands
```bash
git init                    # Initialize repository
git clone <url>            # Clone remote repository
git status                 # Check status
git add <file>             # Stage file
git commit -m "message"    # Commit changes
git push origin main       # Push to remote
git pull origin main       # Pull from remote
```

### Branch Commands
```bash
git branch                 # List branches
git checkout <branch>      # Switch branch
git checkout -b <branch>   # Create and switch
git merge <branch>         # Merge branch
git branch -d <branch>     # Delete branch
```

### History Commands
```bash
git log                    # View commit history
git log --oneline         # Compact history
git diff                  # See changes
git show <commit>         # Show commit details
```

## Troubleshooting Common Issues

### Issue 1: "Repository not found"
```bash
# Check remote URL
git remote -v

# Update remote URL
git remote set-url origin https://github.com/username/repo.git
```

### Issue 2: Merge Conflicts
```bash
# When you see conflict markers, edit the file
# Remove <<<<<<< HEAD, =======, >>>>>>> markers
# Keep the content you want

# Then add and commit
git add conflicted-file.txt
git commit -m "Resolve merge conflict"
```

### Issue 3: Accidental Commits
```bash
# Undo last commit but keep changes
git reset --soft HEAD~1

# Undo last commit and discard changes (careful!)
git reset --hard HEAD~1
```

## Next Steps

Once you're comfortable with these basics:

1. **Learn about Git workflows** (GitFlow, GitHub Flow)
2. **Explore advanced features** (rebasing, cherry-picking, stashing)
3. **Use Git hosting platforms** (GitHub, GitLab, Bitbucket)
4. **Integrate with development tools** (IDEs, continuous integration)

## Summary

Git is an essential tool for any developer. The key concepts to remember:

- **Repository**: Your project folder with Git tracking
- **Commit**: A snapshot of your project at a point in time
- **Branch**: A parallel version of your project
- **Remote**: A copy of your repository on a server

The basic workflow is:  
1. Make changes to files  
2. Stage changes with `git add`  
3. Commit changes with `git commit`  
4. Push to remote with `git push`  

Practice these commands regularly, and you'll become proficient with Git in no time!