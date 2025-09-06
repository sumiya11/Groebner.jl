#!/bin/bash
julia --project=docs -ie 'using Groebner, LiveServer; servedocs()'
