@echo off
for %f in (".\*.eps") do (epstopdf "%f")