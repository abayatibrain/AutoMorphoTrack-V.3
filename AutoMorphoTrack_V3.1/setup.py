
from setuptools import setup, find_packages

setup(
    name='AutoMorphoTrack',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'opencv-python',
        'scikit-image'
    ],
    author='Your Name',
    description='Automated Mitochondria and Lysosome Morphology and Tracking Package',
)
