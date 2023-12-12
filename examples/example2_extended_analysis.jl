# Work in progress

J = V/(nD) # n in rev/sec
n = rpm*(1/60)
J = 0.0
vc = J*n*diameter

scatter!(
    [rpm], [torque], 
    label="RPM = $rpm",
    xlabel="RPM", ylabel="Torque [N]"
    )