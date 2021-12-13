from .generic import ensure_packages_are_installed, update_json_file_kwargs


def setup_colab_runtime_environment(environment_variables_config):
    """
    This function setups up the runtime environment for running the exabyte-api-examples
    notebooks within Google Colaboratory

    Args:
        environment_variables_config (dict): contains key value pairs needed to set up a
            certain notebook (or frontend) runtime.
            Ex) environment_variables_config['ACCOUNT_ID'] = ACCOUNT_ID
                environment_variables_config[notebook_environment] = "Jupyter"

    Returns
        None
    """
    notebook_environment = environment_variables_config['notebook_environment']
    ensure_packages_are_installed(notebook_environment)
    kwargs = {key: environment_variables_config[key] for key in ["ACCOUNT_ID", "AUTH_TOKEN", "MATERIALS_PROJECT_API_KEY", "ORGANIZATION_ID"]}
    update_json_file_kwargs(**kwargs)
