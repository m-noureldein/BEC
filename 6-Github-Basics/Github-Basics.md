# Basic Github Commands to Commit {certain} folder

    #open terminal window and go to the project directory
    cd {path_to_folder}

# Initiate git and push to GitHub

    git config --global user.name "user_name"
    git config --global user.email "user_email"
    git config --global color.ui true
    git init
    git status

# Once you update the folder, add it to git and update the github repository as follows:

    git add *
    git add *
    git commit -a -m "initial commit"
    git status # should state nothing to commit
    git remote add origin  https://{user_name}:{TOKEN}@github.com/{user_name}/{repo}.git 
    git push origin master

# If you get authentication error, get a token by going to GitHub accounts –\> Developer settings –\> personal access token “check all the boxes” –\> generate new token and copy and paste it in the following command:

    git remote add origin https://{user_name}:{TOKEN}@github.com/{user_name}/{repo}.git 
    # or git remote set_url origin https://{user_name}:{TOKEN}@github.com/{user_name}/{repo}.git 
    git remote -v
    git push origin master

# In case none of the above worked, do resetting as follows:

    git fetch --all
    git reset --hard origin/master

# In case you want to create ssh key when you run it on a new computer

    ssh -T -ai ~/.ssh/id_ed25519 git@github.com
    ssh -v -o "IdentitiesOnly=yes" -i ~/.ssh/id_ed25519 git@github.com
    ssh-keygen -t ed25519 -C "user_email" # the email used to create github repository
