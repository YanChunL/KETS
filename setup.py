from setuptools import setup, find_packages

setup(
    name='KETA',
    version='1.0.0',
    author='liyanchun',
    author_email='2812028214@qq.com',
    description='A genomic scaffolding assement tools including Completeness and Accuracy',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/YanChunL/KETA.git',
    packages=find_packages(),
    scripts=['KETA.py'],
    python_requires='>=3.7',
    install_requires=[
        'pandas',
        # Add any other dependencies here
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
