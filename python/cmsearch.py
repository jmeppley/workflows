"""
Suppport methods for running cmsearch
"""
import os
import re


def get_cm_model_version(models_dir):
    """
    Gets the version of the cm and hmm models used from a "VERSION"
    file in the models folder.

    Returns "Unknown" if the version file is missing.
    """

    # Look for a VERSION file
    version_file = os.path.join(models_dir, "VERSION")
    if os.path.exists(version_file):
        with open(version_file) as VF:
            return next(VF).strip()

    # Fall back to README file
    readme_file = os.path.join(models_dir, "README")
    if os.path.exists(readme_file):
        with open(readme_file) as RMF:
            for line in RMF:
                if re.search(r"(Release\s*[0-9.]+)", line):
                    return "RFAM " + line.strip()

    # Fallback to unknown
    return "Unknown"
