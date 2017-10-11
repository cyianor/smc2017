# Workshop and intensive course: Sequential Monte Carlo methods 

## Background

This [workshop](https://www.it.uu.se/conferences/smc2017) on Sequential Monte Carlo (SMC) methods and applications was held from the 30th August to the 1st September 2017 at Uppsala University. Preceding it, a preparatory intensive course concerning the basics of SMC was given from 24th August to 29th August 2017.

The purpose of this repository is to collect course material as well as solutions to the optional exercises which were presented during the course. Since no official solutions exist, this repository is also meant to be a platform for contribution of solutions from other course participants, possibly in a different language than Python (which is what I used). 

The issues section could be used for discussions and questions while the wiki can be used to e.g. summarize discussions which were solved in the issues section. The idea is that this repository should become a valuable resource to course participants who are looking for solutions, as well as a comprehensive source of introductory material on SMC methods.

## Repository structure

Currently, the repository is structured like this:

    smc2017
     |
     |- additional                          Material created during the SMC sessions at FCC
     | 
     |- course_material                     This is where all the distributed material from the course is stored.
     |   |
     |   |- exercise_sheets                 Contains all exercise sheets.
     |   |
     |   |- handwritten_notes               Contains handwritten notes from the course.
     |   |
     |   |- slides                          Contains all slide decks from the course.
     |   |
     |   |- SMC2017.pdf                     The 2017 draft of the SMC lecture notes that Thomas and Fredrik are 
     |                                      working on.
     | 
     |- solutions                           This is where all solutions are stored.
         |
         |- exercises_on_paper              Theoretical exercises solved on paper.
         |
         |- code                            Code for the solution of the applied exercises.

## Guidelines for contribution

Everybody is welcome to clone the repository and send pull requests with additional solutions and notes. 

When you contribute with your own code, it would be great if you could follow the structure below to keep
the repository tidy.

    |- code
        |- Python
        |   |- fheld
        |   |   |- exI.ipynb
        |   |   |- exII.ipynb
        |   |- participant2
        |       |- exI_1.py
        |       |- exI_2.py
        |- Matlab
            |- participant3
                |- exI.m
                |- exII_1and2.m
                |- exII_3and4.m
                |- ouput
                    |- description.txt
                    |- plot1.png
                    |- plot2.png

In summary:
-  Locate the programming language you were using or create a folder for it, in case it does not exist yet
-  Create a folder with some kind of name that separates your solutions from solutions others have written
   in the same language, e.g. `fheld` or `freli`, to avoid confusion about which solutions belong together.
-  File names should follow the pattern
   
        `ex` + exercise sheet roman numeral + file extension (e.g. `.py`)
   
   if the whole exercise sheet is in solved in one file and
   
        `ex` + exercise sheet roman numeral + `_` + exercise numbers (possibly separated by `and`) + file_extension
   
   if only some exercises of an exercise sheet were solved in that file.
-  A folder for output, e.g. plots, can be created in your folder as well. You can either describe in
   your code what output you will get or add an additional file explaining the output.

If Git really isn't your thing (it should though, I feel it is becoming ever more omnipresent), then you can also email me your solutions to felix[dot]held[at]fcc[dot]chalmers[dot]se and I will upload them here.

## Thanks

All the best,

Felix Held (felix[dot]held[at]fcc[dot]chalmers[dot]se)<br>
Fraunhofer-Chalmers Centre<br>
Gothenburg, Sweden

and

Fredrik Lindsten (fredrik[dot]lindsten[at]it[dot]uu[dot]se)<br>
Uppsala University<br>
Uppsala, Sweden
