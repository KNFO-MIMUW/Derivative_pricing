# Derivative_pricing
The aim of the project is to build an application, which allows user to price and determine the sensitivities of different financial derivatives. We plan to use different stochastic, local volatility models and numerical methods to solve our problem. Project will be written in statistical package R and in order to build an application we will use package Shiny. Plan of the project is as follows.
1. We build an application framework in Shiny.
2. We select specific underlying (or set of underlyings) - it is important for us to have access
to volatility surface, because implemented pricing models will be calibrated to volatility
surface.
3. We implement different models that will be calibrated to data choosen in the previous
step. Models to consider:
• Heston model;
• CEV (constant elasticity of variance) model;
• Bates model;
• other models from the class of local and stochastic volatility models.;
In majority of cases we cannot represent price or sensitivities by means of explicit formulas,
hence there is a need to use numerical methods to PDE’s, SDE’s and Monte
Carlo methods (we will try to implement not only basic MC but also advanced variance
reduction methods).
4. Visualize results in the application
