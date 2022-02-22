---
title: S3 Classes
parent: Shared Files
grand_parent: Stage 2 Apps
has_children: false
nav_order: 1
---

## S3 classes

Files in **global/classes** define the object oriented programming (OOP)
components of Stage 2 apps, as R S3 classes.

By convention, a class is defined by scripts named:
    
- **myClass_constructor.R** - for object instantiation
- **myClass_methods.R** - where generic S3 methods are found
- **myClass_utilities.R** - functions called by the above scripts

Constructor functions follow suggested R best practices in being
called 'new_myClass' by convention. Note the case pattern.

Objects defined by classes are never used to access or create
the UI, only to streamline code development by encapsulating 
methods and other common logic. 

### genericMethods.R

If your classes implement non-standard generic S3 methods, 
you must define them in 'genericMethods.R'. 
See comments within for details.

For example, your S3 class can implement the standard generic 
<code>print</code> without worries, but a custom
generic <code>myMethod</code> must be defined in genericMethods.R.
