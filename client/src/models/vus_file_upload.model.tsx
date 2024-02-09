export interface IUnprocessedVus {
  locus: string;
  type: string;
  genotype: string;
  refAllele: string;
  observedAllele: string;
}

export interface IVusGene {
  index: number;
  vus: IUnprocessedVus;
  genes: string[];
}

export interface IVusGeneSelected {
  index: number;
  gene: string;
}

export interface ISamplePhenotypeSelected {
  sampleId: number;
  pheontypeName: string;
  ontologyId: string;
}
