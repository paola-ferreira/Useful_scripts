#### # Getting started in bioinformatics can feel overwhelming. There are countless commands to learn, tools to install, and workflows to understand. But don’t worry. This guide will walk you through the essential things you need before you begin analyzing your data, helping you build confidence step by step.

### ## 1. Getting Access to a Cluster ##

The first step in any bioinformatics is gaining access to a computing cluster. Each institution manages access differently, so the process can vary.  
At Aarhus University (Denmark), for example, access is provided through [GenomeDK](https://genome.au.dk).  
Once your access request has been approved, you can log in and begin working on the cluster environment.

### ## 2. Working With an SSH Client ##

To access a computing cluster, you need an SSH client (a tool that lets you open a secure terminal connection to the cluster).  
If you are using Windows, a highly recommended option is [MobaXterm](https://mobaxterm.mobatek.net).  
I have several reasons why I like this client, but among them, MobaXterm allows you to:  
- Save your login credentials  
- Customize the terminal environment  
- Open multiple SSH sessions  
- Browse remote files through a graphical interface  
- Easily upload and download files using drag-and-drop    
This makes interacting with the cluster much easier, especially if you prefer not to transfer files manually using command-line tools.
If you are on macOS and Linux, the built-in Terminal works well for SSH. However, it does not allow you to interact directly with the files (drag-and-drop).
I am sure that there are SSH clients available for macOS and Linux, but I haven't tested them myself.

### ## 3. Accessing the cluster ##  
To access the cluster, open your SSH client (terminal) and type:  

```ssh user_name@login.genome.au.dk```  

PS: Change to your real username and cluster configurations.  

### ## 4. Creating a Project and a shared folder ##
In genomeDK, you can open a project and share it with your collaborators. It is always a good idea since it is an easier way to share files and ask for help.  
In genomeDK, you can request a project:  

```gdk-project-request -g project__name```

Once your request is approved, you can share with others:  

```gdk-project-promote-user -g <project name> -u username```
