========================
Installation & Updating
========================

Installation
============

First, obtain Python3 (tested on V3.9). If you haven't used python before, we recomend installing it through anaconda.
 `anaconda3 <https://www.anaconda.com/products/individual>`_.

RelaxiFI can be installed using pip in one line. If you are using a terminal, enter:

.. code-block:: python

   pip install RelaxiFI

If you are using Jupyter Notebooks or Jupyter Lab, you can also install it by entering the following code into a notebook cell (note the !):

.. code-block:: python

   !pip install RelaxiFI

You then need to import RelaxiFI into the script you are running code in. In all the examples, we import RelaxiFI as relax.:

.. code-block:: python

   import RelaxiFI as relax

This means any time you want to call a function from RelaxiFI, you do relax.function_name.



Updating
========

To upgrade to the most recent version of RelaxiFI, type the following into terminal:

.. code-block:: python

   pip install RelaxiFI --upgrade

Or in your Jupyter environment:

.. code-block:: python

   !pip install RelaxiFI --upgrade


For maximum reproducability, you should state which version of RelaxiFI you are using. If you have imported RelaxiFI as relax, you can find this using:

.. code-block:: python

    relax.__version__