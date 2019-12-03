def MSR_fuel(nacl,enrichment):
    """Function to give a rough atomic density vector for NaCl-UCl3 salts.
    Parameters
    ----------
    nacl : float
        NaCl percentage in NaCl-UCl3 mixture. Not really mol%.
    enrichment : float
        U235 enrichment in w%
    """
        
    nacl=nacl/100
    enrichment=enrichment/100
    na=nacl/2
    cl=nacl/2+3*((1-nacl)/4)
    u235=((1-nacl)/4)*enrichment
    u238=((1-nacl)/4)*(1-enrichment)
    print(na)
    print(cl)
    print(u235)
    print(u238)
    
