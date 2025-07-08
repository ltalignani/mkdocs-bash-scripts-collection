# Advanced Git: Workflows & Commands

## Introduction

If you have learned [Git Version Control Guide basics](Git.md), it's time to explore advanced concepts that professional developers use daily. This page covers sophisticated Git workflows, powerful commands that can save you time and effort. There is another part that treat [how Git integrates Continuous Integration/Continuous Deployment (CI/CD) systems](ci-cd.md).

These advanced techniques will help you work more efficiently in team environments and handle complex development scenarios.

<h2 class="no-toc">Table of Content</h2>

[TOC]

## Git Workflows

### What is a Git Workflow?

A Git workflow is a recipe or recommendation for how to use Git to accomplish work in a consistent and productive manner. Different teams use different workflows based on their project size, team structure, and deployment strategy.

## GitFlow Workflow

### Overview

GitFlow is a branching model designed around project releases. It's ideal for projects with scheduled releases and provides a robust framework for managing larger projects.

### Branch Types in GitFlow

```
main branch:     A---B---C---D---E---F
                  \                 /
                   \               /
release branch:     G---H----I----J
                     \           /
develop branch:       K----L----M---N---O
                       \       /     \
feature branches:       P-----Q       R---S
                                       \
hotfix branch:                          T---U
```

**Branch descriptions:**

- **main**: *Production-ready* code
- **develop**: Integration branch for features
- **feature**: Individual feature development
- **release**: Preparation for production release
- **hotfix**: Quick fixes for production issues

### GitFlow Commands

First, install GitFlow:
```bash
# On Ubuntu/Debian
sudo apt install git-flow

# On macOS
brew install git-flow-avh
```

#### Initialize GitFlow

```bash
# Initialize GitFlow in your repository
git flow init

# Follow the prompts (usually accept defaults):
# - Production branch: main
# - Development branch: develop
# - Feature branch prefix: feature/
# - Release branch prefix: release/
# - Hotfix branch prefix: hotfix/
```

#### Working with Features

```bash
# Start a new feature
git flow feature start user-authentication

# This creates and switches to: feature/user-authentication
# Equivalent to:
# git checkout develop
# git checkout -b feature/user-authentication

# Work on your feature (add, commit, etc.)
echo "Login form" > login.html
git add login.html
git commit -m "Add login form"

# Finish the feature
git flow feature finish user-authentication

# This merges feature into develop and deletes the feature branch
# Equivalent to:
# git checkout develop
# git merge feature/user-authentication
# git branch -d feature/user-authentication
```

#### Working with Releases

```bash
# Start a release
git flow release start 1.0.0

# This creates release/1.0.0 branch from develop
# Make final adjustments, bug fixes, update version numbers
echo "Version 1.0.0" > VERSION
git add VERSION
git commit -m "Bump version to 1.0.0"

# Finish the release
git flow release finish 1.0.0

# This:
# - Merges release into main
# - Tags the release
# - Merges release back into develop
# - Deletes the release branch
```

#### Working with Hotfixes

```bash
# Start a hotfix (from main branch)
git flow hotfix start critical-bug

# Fix the bug
echo "Bug fix code" > bugfix.py
git add bugfix.py
git commit -m "Fix critical security vulnerability"

# Finish the hotfix
git flow hotfix finish critical-bug

# This:
# - Merges hotfix into main
# - Tags the hotfix
# - Merges hotfix into develop
# - Deletes the hotfix branch
```

### GitFlow Advantages and Disadvantages

**Advantages:**  
- Clear separation of concerns  
- Suitable for scheduled releases  
- Robust for large teams  
- Well-defined process for hotfixes  

**Disadvantages:**  
- Complex for simple projects  
- Can slow down development  
- Requires discipline from team members  

## GitHub Flow

### Overview

GitHub Flow is a simpler alternative to GitFlow, designed for continuous deployment. It's perfect for web applications and projects that deploy frequently.

### GitHub Flow Process

```
main branch:     A---B---D---F-----H
                  \     /     \   /
                   \   /       \ /
feature branch:      C          E
                                 \
another feature:                  G
```

**The process:**

1. **Create a branch** from main
2. **Add commits** to your branch
3. **Open a Pull Request** (PR)
4. **Discuss and review** your code
5. **Deploy** for testing (optional)
6. **Merge** into main

### GitHub Flow Commands

```bash
# 1. Always start from main
git checkout main
git pull origin main

# 2. Create a feature branch
git checkout -b feature/add-search-functionality

# 3. Make changes and commit
echo "Search functionality" > search.js
git add search.js
git commit -m "Add search functionality"

# 4. Push the branch
git push origin feature/add-search-functionality

# 5. Create Pull Request (done on GitHub web interface)
# 6. After review and approval, merge via GitHub
# 7. Clean up locally
git checkout main
git pull origin main
git branch -d feature/add-search-functionality
```

### GitHub Flow Advantages and Disadvantages

**Advantages:**  
- Simple and easy to understand  
- Perfect for continuous deployment  
- Fast development cycle  
- Less overhead than GitFlow  

**Disadvantages:**  
- Less structured for complex projects  
- Requires good testing practices  
- Not suitable for scheduled releases. 

## Advanced Git Commands

### Stashing: Temporarily Save Changes

Stashing allows you to save your current work without committing, useful when you need to quickly switch branches or pull updates.

```bash
# Save current changes to stash
git stash

# Save with a descriptive message
git stash save "Work in progress on login feature"

# List all stashes
git stash list
# Output:
# stash@{0}: WIP on main: 5c3d1a2 Add login form
# stash@{1}: On feature-branch: 8f4e5d6 Work in progress on login feature

# Apply most recent stash
git stash apply

# Apply specific stash
git stash apply stash@{1}

# Apply stash and remove it from stash list
git stash pop

# Show stash contents
git stash show
git stash show -p  # Show full diff

# Delete a stash
git stash drop stash@{1}

# Delete all stashes
git stash clear
```

**Practical example:**
```bash
# You're working on a feature
echo "Incomplete feature" > feature.txt
git add feature.txt

# Emergency: need to fix a bug on main branch
git stash save "Incomplete feature work"

# Switch to main and fix bug
git checkout main
echo "Bug fix" > bugfix.txt
git add bugfix.txt
git commit -m "Fix critical bug"

# Return to feature work
git checkout feature-branch
git stash pop  # Restore your incomplete work
```

### Rebasing: Rewrite History

Rebasing rewrites commit history by moving commits to a new base. It creates a cleaner, linear history.

#### Interactive Rebase

```bash
# Rebase last 3 commits interactively
git rebase -i HEAD~3

# This opens an editor with:
# pick abc123 First commit
# pick def456 Second commit
# pick ghi789 Third commit
```

**Rebase commands:**  
- `pick`: Keep the commit as-is  
- `reword`: Keep commit but edit message  
- `edit`: Keep commit but stop to make changes  
- `squash`: Combine with previous commit  
- `drop`: Remove the commit  
  
**Example - Squashing commits:**
```bash
# Change to:
# pick abc123 First commit
# squash def456 Second commit
# squash ghi789 Third commit

# This combines all three commits into one
```

#### Rebase onto Another Branch

```bash
# Rebase current branch onto main
git rebase main

# Rebase feature branch onto develop
git checkout feature-branch
git rebase develop
```

**Visual representation:**
```
Before rebase:
main:     A---B---C
               \
feature:        D---E---F

After rebase:
main:     A---B---C
                   \
feature:            D'---E'---F'
```

#### Handling Rebase Conflicts

```bash
# If conflicts occur during rebase
git rebase main

# Edit conflicted files, then:
git add conflicted-file.txt
git rebase --continue

# To abort rebase:
git rebase --abort
```

### Cherry-picking: Select Specific Commits

Cherry-picking allows you to apply specific commits from one branch to another.

```bash
# Apply commit abc123 to current branch
git cherry-pick abc123

# Apply multiple commits
git cherry-pick abc123 def456 ghi789

# Apply a range of commits
git cherry-pick abc123..ghi789

# Cherry-pick without committing (stage changes only)
git cherry-pick --no-commit abc123
```

**Practical example:**
```bash
# You have a bug fix in feature branch that's needed in main
git log --oneline feature-branch
# abc123 Fix login bug
# def456 Add new feature (not ready)
# ghi789 Update documentation

# Cherry-pick just the bug fix
git checkout main
git cherry-pick abc123  # Only applies the bug fix
```

**Visual representation:**
```
Before cherry-pick:
main:     A---B---C
               \
feature:        D---E---F (E contains the fix we want)

After cherry-pick:
main:     A---B---C---E'
               \
feature:        D---E---F
```

### Advanced Git Commands Summary

```bash
# Stashing
git stash                    # Save current work
git stash pop               # Restore and remove from stash
git stash list              # Show all stashes

# Rebasing
git rebase -i HEAD~3        # Interactive rebase last 3 commits
git rebase main             # Rebase current branch onto main
git rebase --continue       # Continue after resolving conflicts

# Cherry-picking
git cherry-pick abc123      # Apply specific commit
git cherry-pick abc123..def456  # Apply range of commits
```

## Git Workflow Best Practices

1. **Choose the right workflow** for your team size and deployment frequency
2. **Use descriptive branch names** (`feature/user-auth`, `hotfix/login-bug`)
3. **Keep commits small and focused**
4. **Write clear commit messages**
5. **Use pull requests for code review**

## Advanced Commands Best Practices

1. **Never rebase public branches** (shared with others)
2. **Test after rebasing** before pushing
3. **Use stash for quick context switching**
4. **Cherry-pick sparingly** (can create duplicate commits)
5. **Use interactive rebase** to clean up history before merging


## Conclusion

Advanced Git workflows are essential for professional software development. GitFlow provides structure for complex projects, while GitHub Flow offers simplicity for continuous deployment. Advanced commands like rebasing, cherry-picking, and stashing give you powerful tools for managing your codebase.

Remember: start simple and gradually introduce more advanced practices as your team and project grow. The key is consistency and automation to reduce manual errors and improve code quality.