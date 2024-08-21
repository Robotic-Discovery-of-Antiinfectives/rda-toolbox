Examples
=============

Installation/Usage:
*******************

::

    pip install rda-toolbox

or

::

    pip install "git+https://github.com/Robotic-Discovery-of-Antiinfectives/rda-toolbox.git"


.. code-block:: Python
    :linenos:

    #!/usr/bin/env python3

    import rda_toolbox as rda
    import glob

    rda.readerfiles_rawdf(glob.glob("path/to/raw/readerfiles/*"))
