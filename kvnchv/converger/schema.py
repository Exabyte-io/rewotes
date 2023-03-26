"""Input schema for Manager convergence workflow."""

solver_schema = {
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
        },
    },
    "required": ["name"]
}

input_schema = {
    "type": "object",
    "properties": {
        "solver_input": solver_schema,
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
        },
        "parameter_space": {
            "description": "Parameteric space to evaluate convergence on.",
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {
                        "description": "Parameter name",
                        "type": "string"
                    },
                    "start": {
                        "description": "Initial parameter value.",
                        "type": "number"
                    },
                    "max": {
                        "description": "Final parameter value.",
                        "type": "number"
                    },
                    "delta": {
                        "description": "Step size between parameters"
                    }
                },
                "required": ["name", "start", "max"]
            }
        }
    },
    "required": ["solver_input", "input_path", "tol", "target", "parameter_space"]
}
