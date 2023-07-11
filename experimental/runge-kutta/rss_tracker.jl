using Printf

function setup_memuse_tracker()
    tracker = Ref(0)
    function mem_use(tracker)
        finalizer(mem_use, tracker)
        out = Core.CoreSTDOUT()
        Core.write(
            Core.CoreSTDOUT(),
            @sprintf "GC live: %9.3f MiB,    " Base.gc_live_bytes() / 2^20
        )
        Core.write(Core.CoreSTDOUT(), @sprintf "Max. RSS: %9.3f MiB\n" Sys.maxrss() / 2^20)
        nothing
    end

    finalizer(mem_use, tracker)
    nothing
end
