"""Input schema for Managerconvergence workflow."""
input_schema = {
    "type": "object",
    "properties": {
        "solver_input": {
            "description": "Solver parameters",
            "type": "object",
            "properties": {
                "name": {
                    "description": "Name of solver type.",
                    "type": "string"
                },
                "solver_path": {
                    "description": "Explicit path to solver.",
                    "type": "string"
                }
            },
            "required": ["name"]
        },
        "input_path": {
            "description": "Path containing required files for solver.",
            "type": "string"
        },
        "tol": {
            "description": "Minimum tolerance for target quantity.",
            "type": "number"
        },
        "target": {
            "description": "Output quantity to converge.",
            "type": "string",
            "enum": [
                "total_energy"
            ]
        },
    },
    "required": ["solver_input", "input_path", "tol", "target"]
}
