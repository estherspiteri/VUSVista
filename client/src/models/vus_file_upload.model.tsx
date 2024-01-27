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
  genes: string;
}
