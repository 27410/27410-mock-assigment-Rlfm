# 27410 - Mock exam

In this mock exam you will solve a number of tasks related to topics covered in class. It is supposed to give you an idea of how the final exam is going to be. You should be able to solve it in <3 hours (please track how much time you spent on it).

In `assignment.py`:

1. Based on the model, what is the theoretical maximum yield of acetate in units of mmol-acetate/mmol-glucose?
2. Based on the model, what is the theoretical maximum yield of acetate in units of cmol-acetate/cmol-glucose?
3. Based on the model's stoichiometry alone, how many reaction fluxes need to be measured in order to make the system determined and solvable?
4. How much is the (optimal) growth rate reduced if fumerase (FUM in the model) is overexpressed to have a 2-fold higher flux in comparison to its flux at maximum growth rate?
5. What genes are essential under acetate conditions but not glucose conditions?

Hints:
* `assignment.py` contains hints and recommendations both in the beginning and in the comments to each task.
* Remember to undo modifications to the model before continuing with the next task (either make a copy of the model for each task or use the `with model: ...` statement as shown in the exercise).
* If not mentioned otherwise, use the default conditions set in the model in the beginning of the file.

Jupyter classroom backup options:
* Option 1: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/27410/course-materials/master?urlpath=lab) You can work on the exam/assignment by clicking on the binder badge to the left. You will have to use git to clone the exam/assignment repository before you can get startet (also you will have to configure your email and name the first time you commit). **WARNING no work is saved on the Binder Jupyterhub so back up your work regularly with git**
* Option 2: Follow the installation instructions [here](https://github.com/27410/course-materials/blob/master/INSTALLATION.md) to install git, Python and cobrapy on your local computer.

*The Rules*:
* Every tasks needs to be solved programmatically going from the provided data to the final result, do not use manual steps in between.
* The ideal solution contains code only, is self-explanatory and provides sufficient context for the teacher to follow the solution process.
* Do not replace the placeholders _ with numbers that you took directly out of the test cases, manually calculated or guessed. Ideally placeholders _ should be replaced with the final calculation step or a variable that contains the result of the final calculation step.
* Use human readable variable/function names (avoid cryptic abbreviations and one letter variables). Consider the recommendations laid out in https://www.python.org/dev/peps/pep-0008/
