# Install stable version through PyPi
There are a few [dependencies](http://hichipper.readthedocs.io/en/latest/content/Dependencies.html)
needed to get **hichipper** to run. All are 
very common bioinformatics tools / languages and should be readily available in
most systems. However, **note that the current implementation of hichipper is not supported
on Windows platforms**. 

Depending on your python environment, we generally recommend using a virtual environment
to keep python dependencies tidy. An example of installing **hichipper** inside a new
python virtual environment called `venv` using the following sequence of commands--

```
virtualenv -p /usr/bin/python2.7 venv
source venv/bin/activate
pip install hichipper
hichipper --version
```

# Install via GitHub

Though **not recommended**, a bleeding-edge (development) version can be installed
directly from Git. Again using a virtual environment--

```
virtualenv -p /usr/bin/python2.7 venv
source venv/bin/active
pip3 install git+ssh://git@github.com/aryeelab/search/tree/master/hichipper
```

While installing **hichipper** is obviously a great first step, make sure that all of the 
[dependencies](http://hichipper.readthedocs.io/en/latest/content/Dependencies.html) are met. 
Check out the next page for more detail. 
<br><br>