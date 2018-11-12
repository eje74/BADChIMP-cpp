Notes on the Makefile:  
* Note on file name: `Makefile.terminal` is used since my IDE already creates an `Makefile`.  To use this with the make command write  
`make -f Makefile.terminal`
* The `TARGET` is the main executable of the project.
Type make and this is what gets built.
* I use the fewest number of compiler `CFLAGS` when developing as possible,
optimization happens later.
* The `SOURCES` list is dynamic, I don’t want to manually have to maintain this list as I program. Anything in the `src` folder will be included in the compile as long as it has a `SRCEXT` extension.
* The `OBJECTS` list is also dynamic and uses a Makefile trick to build the
list based on available sources.
* The `LIB` add local(?) libraries.
* The `INC` ensures all headers in the include folder are accessible.
* He likes to see the commands that run, hence the multitude of `@echo`'s.
* Manually add snippets and test builds as a new Makefile target.
* `.PHONY`: A phony target is one that is not really the name of a file; rather it is just a name for a recipe to be executed when you make an explicit request. There are two reasons to use a phony target: to avoid a conflict with a file of the same name, and to improve performance.

Run
```
mkdir -p bin build doc include lib snippets src test
```
to build the directory structure.

Based on a blog-post by *Hilton Lipschitz*
