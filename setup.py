from setuptools import setup, find_packages

setup(name='VHTS',
        version='0.1',
        packages=['vhts'],
#        packages=find_packages(),
        url='https://github.com/gicsaw/VHTS',
        license='MIT LICENSE',
        author='Seung Hwan Hong',
        author_email='gicsaw0@gmail.com',
        description='',
        scripts=['bin/master_dock.py',
                 'pydock_run.py',
                 'sub_dock.py',
                 'vhts_check_restart.py'
                ]
)


