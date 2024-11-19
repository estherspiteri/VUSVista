export interface IClinvarClassificationUpdate {
  classification: string;
  reviewStatus: string;
  lastEval: string;
}

export interface IClinvarUpdate {
  dateChecked: string;
  update?: IClinvarClassificationUpdate | null;
}
